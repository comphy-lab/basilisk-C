#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "tag.h"
#include "navier-stokes/perfs.h"
#include "timedependent_bcs.h"
#include "view.h"

// Geometry
double L = 6.5e-3;
double W = 4e-3;
double ys = 0.95e-4;         
double ye = 2e-4;

double length = 0.5e-3;

// Prefilmer geometry
double xsP = 0.;
double xeP = 1e-3;
double w = 200e-6;           

// Grid parameters
double uemax = 0.5;
int maxlevel = 9;
double eps = 1e-6; 

double SIGMA = 0.0275;
double u0 = 0.5;
double uAir = 50;
//BCs
InletData inletTop = {0};

#define IN_LIQ (yy >= ys && yy <= ye)
#define IN_WALL (yy < ys && yy >= -w/2)
#define IN_AIR_UP (yy > ye)
#define IN_AIR_DOWN (yy < -w/2.)

//BCs
scalar f0[];

static double inlet_ux (double yy, double zz, double tt)
{
  if IN_WALL
    return 0;
  else if IN_AIR_UP
    return sample_inlet_component (&inletTop, inletTop.Ux, tt, yy, zz);
  else if IN_AIR_DOWN 
    return sample_inlet_component (&inletTop, inletTop.Ux, tt, -yy+1e-4, zz);
  else if IN_LIQ
    return u0;
  else
    return 0.;
}

static double inlet_uy (double yy, double zz, double tt)
{
  if IN_WALL
    return 0;
  else if IN_AIR_UP
    return sample_inlet_component (&inletTop, inletTop.Uy, tt, yy, zz);
  else if IN_AIR_DOWN 
    return -sample_inlet_component (&inletTop, inletTop.Uy, tt, -yy + 1e-4, zz);
  else if IN_LIQ
    return 0.;
  else
    return 0.;
}

static double inlet_uz (double yy, double zz, double tt)
{
  if IN_WALL
    return 0;
  else if IN_AIR_UP
    return sample_inlet_component (&inletTop, inletTop.Uz, tt, yy, zz);
  else if IN_AIR_DOWN 
    return sample_inlet_component (&inletTop, inletTop.Uz, tt, -yy + 1e-4, zz);
  else if IN_LIQ
    return 0.;
  else
    return 0.;
}

static double inlet_f (double yy, double zz, double tt)
{
  return (yy >= ys && yy <= ye) ? 1. : 0.;
}


// ========= Left =========
u.n[left]  = dirichlet(inlet_ux(y,z,t));
u.t[left]  = dirichlet(inlet_uy(y,z,t));
u.r[left]  = dirichlet(inlet_uz(y,z,t));
//uf.n[left] = inlet_ux(y,z,t);

p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet (inlet_f(y, z, t));

// =========    Right    =========//
u.n[right]  = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));
u.t[right]  = neumann(0.);
u.r[right]  = neumann(0.);
uf.n[right] = neumann(0.);

p[right]  = dirichlet(0.);
pf[right]  = dirichlet(0.);
f[right]  = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));

// ===========   Top   =========== //
u.n[top]  = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));
u.t[top]  = neumann(0.);
u.r[top]  = neumann(0.);

p[top]    = dirichlet(0.);
pf[top]   = dirichlet(0.);
f[top]    = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));

// ==========  Bottom   ========== //
u.n[bottom]  = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));
u.t[bottom]  = neumann(0.);
u.r[bottom]  = neumann(0.);

p[bottom]    = dirichlet(0.);
pf[bottom]   = dirichlet(0.);
f[bottom]    = (uf.n[] > 0. ? neumann(0.) : dirichlet(0.));

// =======  Embed  ====== //
u.n[embed] = dirichlet(0.);
u.t[embed] = (fabs(z) > W/2 - eps ? neumann(0.) : dirichlet(0.));
u.r[embed] = (fabs(z) > W/2 - eps ? neumann(0.) : dirichlet(0.));
f[embed] = neumann(0.);
p[embed]  = neumann(0.);
pf[embed]  = neumann(0.);


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

int main()
{
    init_grid(1 << 8);
    origin( 0. , -L/2 , -L/2);
    size(L);

    mu1 = 0.0015631;
    mu2 = 1.8e-5;
    rho1 = 770;
    rho2 = 1.2;

    f.sigma = SIGMA;
    G.y = -9.81;
    TOLERANCE = 2e-4;
    run();
}

event init(i=0)
{   
    if (pid() == 0) {
        int ret = system("gunzip -k -f inletU_cut.bin.gz");
    }
    #if _MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    if (load_inlet_data("inletU_cut.bin", &inletTop) != 0)
    exit(1);
    
    CFL = 0.6;
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

// Numerical
event adapt(i+=5)
{
    adapt_wavelet({cs,f,u},(double[]){1e-5,0.01,uemax,uemax,uemax},maxlevel);
    build_geometry();
}

event sponge (i++) {
  double tau = 4e-6;  
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

event snapshot (t = 0; t += 1e-4; t <= 10e-3) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  scalar pid[];
  foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());
  dump (name);
}

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


