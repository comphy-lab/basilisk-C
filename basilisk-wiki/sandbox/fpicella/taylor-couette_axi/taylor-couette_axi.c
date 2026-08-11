/**
# Taylor--Couette flow around a hollow rotating cylinder

We consider a rotating hollow cylindrical shell of outer radius $R_1$
inside a fixed cylinder of radius $R_2$, with fluid on both sides of the
shell. The radius ratio is $\eta=R_1/R_2\simeq0.5$, and

$$
Re_i=\frac{\Omega R_1(R_2-R_1)}{\nu}.
$$

For $\eta=0.5$, circular Couette flow becomes linearly unstable at
$Re_{i,c}\simeq68$ [Esser & Grossmann (1996)](#references). Here
$\Omega=2$ and $\nu=0.02$, so $Re_i\simeq100$ and axisymmetric Taylor
vortices are expected.

![Pressure field, meridional streamlines and velocity vectors during
the development of Taylor vortices.](taylor-couette_axi/pressure-velocity.mp4)

The rotating cylinder is hollow: the embedded solid shell occupies
$R_1-t<r<R_1$, while fluid is also solved for in the inner region
$r<R_1-t$.
*/

#include "grid/quadtree.h"
#include "axi.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "swirl_embed.h"
#include "view.h"
#include "axistream.h"

double endTime = 100.;

// ============================================================
// Geometry
// ============================================================

double R1 = 1. + 1e-5;
double R2 = 2.;

double thickness = 0.1;

#define RIN (R1 - thickness)

// ============================================================
// Fluid properties
// ============================================================

double nu = 0.02;
double Omega = 2.;

// ============================================================
// Resolution
// ============================================================

int LEVEL = 6;

// ============================================================
// Viscosity
// ============================================================

face vector muv[];

// ============================================================
// Axis r = 0
// ============================================================

u.n[bottom] = dirichlet (0.);
u.t[bottom] = neumann (0.);

w[bottom] = dirichlet (0.);

// ============================================================
// Fixed outer cylinder r = R2
// ============================================================

u.n[top] = dirichlet (0.);
u.t[top] = dirichlet (0.);

w[top] = dirichlet (0.);

// ============================================================
// Rotating embedded shell
// ============================================================

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

w[embed] = dirichlet (Omega*y);

// ============================================================
// Embedded geometry
// ============================================================

void shell_geometry()
{
  vertex scalar phi[];

  foreach_vertex()
    phi[] = (y - RIN)*(y - R1);

  boundary ({phi});

  fractions (phi, cs, fs);

  fractions_cleanup (cs, fs);
}

// ============================================================
// Main
// ============================================================

int main()
{
  size (R2);

  origin (0., 0.);

  // Periodic axial direction.
  periodic (right);

  init_grid (1 << LEVEL);

  mu = muv;

  DT = 1.e-1;

  TOLERANCE = 1.e-6;

  display_control (Omega, 0., 10.);

  run();
}

// ============================================================
// Properties
// ============================================================

event properties (i++)
{
  foreach_face()
    muv.x[] = nu*fm.x[];

  boundary ((scalar *) {muv});
}

// ============================================================
// Initialization
// ============================================================

event init (t = 0)
{
  shell_geometry();

  double eps = 1.e-4;

  // One axial perturbation wavelength over the periodic domain.
  double k = 2.*pi/L0;

  foreach() {

    // Base azimuthal flow.
    w[] = Omega*y;

    u.x[] = 0.;
    u.y[] = 0.;

    // Perturb the fluid on both sides of the shell.
    if (cs[] > 0.) {
      u.x[] += eps*cos(k*x);
      u.y[] += eps*sin(k*x);
      w[]   += eps*cos(k*x);
    }
  }

  boundary ((scalar *) {u, w});
}


/**
## Development of Taylor vortices

Pressure is shown in colour, white lines are meridional streamlines and
arrows represent the meridional velocity $(u_z,u_r)$.

*/

event movie_pressure_velocity (t += 0.5)
{
  scalar psi[];

  axistream (u, psi);

  stats spsi = statsf (psi);

  double psimax =
    max (fabs(spsi.min), fabs(spsi.max));

  if (psimax < 1.e-12)
    psimax = 1.e-12;

  /*
   * Pressure range at the current time.
   */
  stats sp = statsf (p);

  double pmin = sp.min;
  double pmax = sp.max;

  view (
    fov = 25,
    tx = -0.5,
    ty = -0.5,
    width = 900,
    height = 900
  );

  clear();

  // Pressure background.
  squares (
    "p",
    linear = true,
    min = pmin,
    max = pmax,
    cbar = true,
    label = "p"
  );

  // Meridional streamlines.
  isoline (
    "psi",
    n = 21,
    min = -psimax,
    max =  psimax,
    lw = 0.5,
    lc = {1., 1., 1.}
  );

  // Meridional velocity.
  vectors (
    "u",
    scale = 0.05
  );

  // Rotating shell.
  draw_vof (
    "cs",
    "fs",
    lw = 3.
  );

  box();

  char label[200];

  sprintf (
    label,
    "t = %.2f",
    t
  );

  draw_string (
    label,
    pos = 1,
    size = 18
  );

  save ("pressure-velocity.mp4");
}


/**
## Pressure on the embedded shell

At the final time, pressure is interpolated directly at the barycentre
of each embedded-boundary fragment. The two sets correspond to the
inner and outer faces of the hollow rotating cylinder.
*/

event pressure_profiles (t = endTime)
{
  FILE * fi = fopen ("pressure-inner.dat", "w");
  FILE * fo = fopen ("pressure-outer.dat", "w");

  fprintf (fi, "# z p\n");
  fprintf (fo, "# z p\n");

  foreach() {

    if (cs[] > 0. && cs[] < 1.) {

      /*
       * Position and normal of the embedded fragment
       * relative to the cell centre.
       */
      coord b, n;

      double area = embed_geometry (point, &b, &n);

      if (area > 0.) {

        /*
         * Physical position of the embedded-fragment
         * barycentre.
         */
        double zb = x + b.x*Delta;
        double rb = y + b.y*Delta;

        /*
         * Pressure interpolated directly at the embedded
         * boundary.
         */
        double pb = embed_interpolate (point, p, b);

        /*
         * Identify which side of the hollow shell this
         * fragment belongs to.
         */
        if (fabs(rb - RIN) < fabs(rb - R1))
          fprintf (fi, "%.12g %.12g\n", zb, pb);
        else
          fprintf (fo, "%.12g %.12g\n", zb, pb);
      }
    }
  }

  fclose (fi);
  fclose (fo);
}


/**
The pressure profiles along both embedded surfaces are shown below.
`smooth unique` also orders the points in the axial direction before
drawing the curves.

~~~gnuplot Pressure along the rotating shell
set xlabel 'z'
set ylabel 'p'
set grid
set key top right

plot \
  'pressure-inner.dat' using 1:2 smooth unique \
    with lines lw 2 title 'inner surface', \
  'pressure-outer.dat' using 1:2 smooth unique \
    with lines lw 2 title 'outer surface'
~~~

## References

* A. Esser and S. Grossmann, *Analytic expression for Taylor--Couette
  stability boundary*, Physics of Fluids **8**, 1814--1819 (1996).

* S. Grossmann, D. Lohse and C. Sun, *High-Reynolds Number
  Taylor--Couette Turbulence*, Annual Review of Fluid Mechanics **48**,
  53--80 (2016).
*/


// ============================================================
// Stop
// ============================================================

event stop (t = endTime)
{
  return 1;
}