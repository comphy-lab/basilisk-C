/**
# Planar prefilmer airblast atomizer

This test case models the interaction between a liquid film and a
high-speed coflowing air stream in a planar prefilmer configuration.

The incompressible two-phase Navier--Stokes equations are solved on an
adaptive octree grid with surface tension and embedded boundaries. The
conserving formulation is used to improve momentum conservation, while
`tag()` is used to identify connected liquid structures for statistical
post-processing.
*/

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "tag.h"
#include "navier-stokes/perfs.h"
#include "view.h"

/**
We define the domain size and a domain width that will be cut on the z direction with embed. 
The prefilmer is represented by a thin solid element located between `xsP`
and `xeP`. Its thickness is set by `w`. This geometry affects the formation
and transport of the liquid film before breakup.
*/

// Geometry
double L = 6.5e-3;
double W = 4e-3;
double ys = 1e-4;         
double ye = 2e-4;

// Prefilmer geometry
double xsP = 0;
double xeP = 1e-3;
double w = 200e-6;           

/**

The maximum grid resolution is controlled by `maxlevel`, while `eps` is
a small tolerance used in the treatment of embedded boundaries. The
parameter `length` defines the initial streamwise extent of the liquid
film.

The reference velocity `uemax` is used in the mesh adaptation criterion.
Surface tension is set by `SIGMA`. The liquid and gas enter the domain
with velocities `u0` and `uAir`, respectively.
*/

double uemax = 0.5;
int maxlevel = 9;
double eps = 1e-6; 
double length = 0.2e-3;

double SIGMA = 0.0275;
double u0 = 0.5;
double uAir = 50;

scalar f0[];

/**
## Boundary Conditions

To simplify the prescription of inlet boundary conditions, we distinguish
three regions at the left boundary:

- `IN_WALL`: the prefilmer wall region,
- `IN_AIR`: the gas inlet region,
- `IN_LIQ`: the liquid inlet region.

These masks are used to prescribe the inlet velocity and phase fraction as functions of the vertical coordinate.
The right, top, and bottom boundaries are treated as open outlets with no recirculation, whereas the lateral walls are modeled as slip walls. Dirichlet conditions are imposed at the left boundary.
*/

#define IN_WALL (y <= w/2 && y >= -w/2)
#define IN_AIR ( y > ye || y < -w/2)
#define IN_LIQ ( y >= ys && y <= ye)

// Left 
u.n[left] = dirichlet ( IN_LIQ ? u0 : IN_AIR ? uAir  : IN_WALL? 0. : 0.);
u.t[left] = dirichlet(0.);
u.r[left] = dirichlet(0.);
p[left]   = neumann(0.);
f[left] = dirichlet (IN_LIQ ? 1. : IN_WALL ? 0. : 0.);

// Right   
u.n[right]  = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));
u.t[right]  = neumann(0.);
u.r[right]  = neumann(0.);
uf.n[right] = neumann(0.);

p[right]  = dirichlet(0.);
pf[right]  = dirichlet(0.);
f[right]  = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));

// Top   
u.n[top]  = neumann(0.);
u.t[top]  = neumann(0.);
u.r[top]  = neumann(0.);

p[top]    = dirichlet(0.);
pf[top]   = dirichlet(0.);
f[top]    = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));

// Bottom 
u.n[bottom]  = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));
u.t[bottom]  = neumann(0.);
u.r[bottom]  = neumann(0.);

p[bottom]    = dirichlet(0.);
pf[bottom]   = dirichlet(0.);
f[bottom]    = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));

// Embed  
u.n[embed] = dirichlet(0.);
u.t[embed] = (fabs(z) > W/2 - eps ? neumann(0.) : dirichlet(0.));
u.r[embed] = (fabs(z) > W/2 - eps ? neumann(0.) : dirichlet(0.));
f[embed] = neumann(0.);
p[embed]  = neumann(0.);
pf[embed]  = neumann(0.);

/**
## Geometry

The computational domain has streamwise length `L` and spanwise width `W`.
The spanwise extent is enforced using embedded boundaries.

The prefilmer is modeled as a thin solid plate extending from `xsP` to
`xeP`, with thickness `w`. The liquid inlet is located between `ys` and
`ye`.
*/

static inline double phi_prefilmer_solid (double x, double y)
{
  return min( min( x - xsP , xeP - x ), w/2 - fabs(y) );
}

static inline double phi_channel (double z)
{
  return W/2 - fabs(z); 
}

static void build_geometry (void)
{
  vertex scalar phi[];

  foreach_vertex() {
    double solid = phi_prefilmer_solid(x, y);
    double chan  = phi_channel(z);           

    phi[] = min(chan, -solid);
  }

  boundary({phi});
  fractions (phi, cs, fs);
  cs.refine = cs.prolongation = fraction_refine;
  fractions_cleanup (cs, fs);
  restriction ({cs, fs});
}

/**
## Main parameters

The liquid and gas properties are chosen to match the experimental
conditions. Surface tension is set through `f.sigma`, gravity acts in
the negative `y` direction, and the pressure solver tolerance is reduced
to improve stability.
*/
int main()
{
    init_grid(1 << 9);
    origin( 0. , -L/2 , -L/2);
    size(L);

    mu1 = 0.0015631;
    mu2 = 1.8e-5;
    rho1 = 770;
    rho2 = 1.2;

    f.sigma = SIGMA;
    G.y = -9.81;
    TOLERANCE = 1e-4;
    run();
}
/**
## Initial conditions

The embedded geometry is constructed first. The liquid phase is then
initialized within the inlet region and limited to a finite streamwise
extent. The axial velocity is finally initialized consistently with the
local phase distribution, with `u0` in the liquid and `uAir` in the gas.
*/
event init(i=0)
{   
    build_geometry();
    fraction(f0, min(min(y - ys, ye - y),min(z + W/2, W/2 - z)));
    f0.refine = f0.prolongation = fraction_refine;
    restriction({f0});

    foreach() 
    {
      // Setup f
      f[] = cs[] * f0[] * (x < length);
      f[] = clamp (f[], 0., 1.);

      if (cs[] > 0.) {

          double velocity_fluid = u0 * f[] + uAir * (1. - f[]);
          u.x[] = velocity_fluid * cs[];
      } 
    }
}

/**
## Mesh adaptation

The mesh is dynamically adapted using wavelet-based error estimates on
the embedded geometry, volume fraction, and velocity field.
*/

event adapt(i++)
{
    adapt_wavelet({cs,f,u},(double[]){1e-5,0.01,uemax,uemax,uemax},maxlevel);
    build_geometry();
}

/**
## Sponge layer

A damping layer is applied near the outlet in order to reduce spurious
reflections and stabilize the outflow. In this region, the velocity field
is progressively relaxed toward a uniform gas stream with axial velocity
`uAir` and zero transverse components.

The damping intensity increases smoothly from `x0` to the outlet over a
characteristic timescale `tau`.
*/

event damp (i++) {
  double tau = 1e-6;  
  double x0  = 0.8*L;  // start sponge 

  foreach() {
    if (x > x0) {
      double s = (x - x0)/(L - x0);   
      s = clamp(s, 0., 1.);

      s = s*s;

      double a = (dt/tau)*s;         
      if (a > 1.) a = 1.;

      u.x[] += a*(uAir - u.x[]);
      u.y[] += a*(0.   - u.y[]);
      u.z[] += a*(0.   - u.z[]);
    }
  }
  boundary ((scalar*){u});
}

/**
## Output

This event identifies the connected liquid structures and computes their
volume and centroid. The largest connected component is assumed to be
the liquid core, and its streamwise extent is extracted as a function of
the spanwise coordinate.

The resulting core profile is written to file, together with the volume
and centroid of all detected liquid structures. 


*/

event ligaments_and_droplets (t += 1e-4)
{
  scalar m[];
  foreach()
    m[] = f[] > 1e-3;

  int n = tag(m);

  if (n > 0) {

    // Arrays storing the volume and centroid of each structure
    double v[n];
    coord b[n];

    for (int j = 0; j < n; j++) {
      v[j] = 0.;
      b[j].x = b[j].y = b[j].z = 0.;
    }

    // Accumulate volume and centroid coordinates
    foreach (serial)
      if (m[] > 0) {
        int j = m[] - 1;
        double vol = dv()*f[];

        v[j] += vol;
        b[j].x += vol*x;
        b[j].y += vol*y;
        b[j].z += vol*z;
      }

#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    // Identify biggest connected structure
    double highest = 0.;
    int core_tag = -1;

    for (int j = 0; j < n; j++) {
      if (v[j] > highest) {
        highest = v[j];
        core_tag = j + 1; 
      }
    }

    // Discretization along z for profile
    int n_bins = 400;
    double dz = L0/n_bins;
    double x_max[n_bins];

    for (int k = 0; k < n_bins; k++)
      x_max[k] = -1e30;

    // Find x max for each bin z of the core
    foreach (serial)
      if (m[] == core_tag) {
        int idx = (z - Z0)/dz;
        if (idx >= 0 && idx < n_bins && x > x_max[idx])
          x_max[idx] = x;
      }

#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, x_max, n_bins, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    // Write core profile
    if (pid() == 0) {
      char name[80];
      sprintf(name, "profilo_core_z-%g.dat", t);
      FILE *fp = fopen(name, "w");

      for (int k = 0; k < n_bins; k++) {
        if (x_max[k] > -1e29) {
          double z_center = Z0 + (k + 0.5)*dz;
          fprintf(fp, "%g %g\n", z_center, x_max[k]);
        }
      }
      fclose(fp);
    }

    // Output droplets statistics 
    for (int j = 0; j < n; j++)
      fprintf(fout, "%d %g %d %g %g %g %g\n",
              i, t, j,
              v[j],
              b[j].x/v[j], b[j].y/v[j], b[j].z/v[j]);

    fflush(fout);
  }
}

/**
Finally, snapshots of the full simulation state are saved every `1e-4` seconds for restart and post-processing purposes.
*/
event snapshot ( t += 1e-4 ) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  scalar pid[];
  foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());
  dump (name);
}

/**
Two interface movies are generated during the simulation. The first one
shows a top view of the liquid structure, while the second one provides
a side view. Both animations are saved every `1e-5` seconds and are used
to monitor the evolution of the interface and the atomization process.
*/
event movie (t += 1e-5) {

  clear();
  view (quat = {-0.433, -0.556, -0.437, 0.559},
        fov = 5, near = 0.01, far = 1000,
        tx = 0.000, ty = -0.489, tz = -11.026,
        width = 1800, height = 1100);
  draw_vof ("f");
  save ("top_view.mp4");
  
  clear();
  view (quat = {-0.002, -0.995, -0.099, -0.017},
        fov = 5, near = 0.01, far = 1000,
        tx = 0.412, ty = 0.077, tz = -11.127,
        width = 1800, height = 1100);
  draw_vof ("f");
  save ("side_view.mp4");
}

// End
event end(t=10e-3)
{
    return 1;
}


