/**
# Vortex shedding behind a sphere at Reynolds $\in [200, 2000]$

This setup is a slight modification of *the* [flow past a sphere
example](http://www.basilisk.fr/src/examples/sphere.c). These
modifications are concerned with generating this movie:

![Flow past a sphere with increasing Reynolds number ($Re$). The
 visualization shows a volumetric rendering of the stagnant fluid
 (red) in the wake, upto the fluid with 90% of the mean velocity
 (blue), with increasing
 transparancy.](https://www.antoonvanhooft.nl/media/sphere.mp4)(width="100%")

This visualization can be compared to the visualization of the $\lambda_2$ citerion

![Vortex detection ($\lambda_2$) isosurface rendering](https://www.antoonvanhooft.nl/media/sphere_l2.mp4)

We solve the Navier--Stokes equations on an adaptive octree and use
embedded boundaries to define the sphere. */

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
// a bwatch movie...
#include "bwatch.h"

/**
We will use the $\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) for vortex detection. */

#include "lambda2.h"

/**
This is the maximum level of refinement i.e. an equivalent maximum
resolution of $1024^3$. */

int maxlevel = 10;

/**
We need a new field to define the viscosity. */

face vector muv[];

/**
The domain size is $16^3$. We move the origin so that the center of
the unit sphere is not too close to boundaries. */

double D = 1. [1];

int main()
{
  init_grid (64);
  size (16.*D);
  origin (- 3.*D, -L0/2., -L0/2.);
  mu = muv;
  run();
}


/**
The viscosity is given by the Reynolds number, the sphere diameter,
the inflow velocity, and fluid fraction.*/

double U0 = 1., Re = 200.;

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Re;
}

/**
The boundary conditions are inflow with unit velocity on the
left-hand-side and outflow on the right-hand-side. */

u.n[left]  = dirichlet(U0);
u.t[left] = dirichlet (0);
u.r[left] = dirichlet (0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The boundary condition is no slip on the embedded boundary. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

/**
$Re$ is increased from 200 to 2000.
 */
event increase_Re (t = 10; t += 1) {
  Re = min (Re*1.015, 2000);
}

event init (t = 0) {

  /**
  We initially refine only in a sphere, slightly larger than the solid
  sphere. The sphere is displace a bit, to make text rendering in the movie more convenient. */

  refine (x*x + sq(y + D) + z*z < sq(1.2*D/2.) && level < maxlevel);

  /**
  We define the unit sphere. */

  solid (cs, fs, x*x + sq(y + D) + z*z - sq(D/2.));

  /**
  We set the initially horizontal velocity to the inflow velocity
  everywhere (outside the sphere). */
  
  foreach()
    u.x[] = cs[] ? U0 : 0.;
}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We use Basilisk view to create the animated isosurface of $\lambda_2$. */

event movies (t += 0.125; t <= 240)
{

  /**
  Here we compute two new fields, $\lambda_2$ and the vorticity
  component in the $y-z$ plane. */
  
  scalar l2[], vyz[];
  foreach()
    vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);
  lambda2 (u, l2);

  view (fov = 11.44, quat = {0.072072,0.245086,0.303106,0.918076},
	tx = -0.307321, ty = 0.22653, bg = {1,1,1},
	width = 802, height = 634);
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.01, color = "vyz", min = -1, max = 1,
	      linear = true, map = cool_warm);
  save ("movie.mp4");
  /**
Here we compute the a flow speed related field, for volumetric rendering with `bwatch`. 
   */
  scalar U[];
  foreach() {
    U[] = 0;
    foreach_dimension()
      U[] += sq(u.x[]);
    U[] = 0.9*U0 - sqrt(U[]);
  }

  static FILE * fp = popen ("ppm2mp4 sphere.mp4", "w");
  double angle = 1.2*cos(2.*pi*t/240.);
  char capt[99];
  sprintf (capt, "Re = %d ", (int)Re);
  watch (fov = 10, O = {1.5*L0*cos(angle), 2*D, 1.5*L0*sin(angle + 1e-6)}, poi = {X0 + L0/2., -D, 1e-3}, nx = 1500, ny = 600);
  // Sketching an exact sphere is cheaper compared to rendering its fov facets
  sphere (D/2.,{1e-6, -D, 0} , mat = {.col = {100, 100, 100}, .R = 0.1});
  // Encapsulating sphere as background
  sphere (4*L0, mat = {.dull = true});
  // Fix me: Is removing these "not exist" files needed?
  system ("rm 3not9exist6.ppm");
  system ("rm 2not8exist5.ppm");
  // Note $Re$ in movie
  sketch_text (capt, ops = "-font Bookman-DemiItalic -rotate 90", res = 2000, 
	       pos = "center", tc = {23, 23, 23}, fs = 40, alpha = 0.01, n = {1, 0, 0});
  // Volumetric rendering
  volume (U, cols = true, sc = 0.25, mval = 0.001, shading = 1);
  store (fp);
  // Fit me: Should image data be handled better (is it even freed?)
  plain();
  // Remainder from Testing
  //store (fopen ("test.ppm", "w"));
}

/**
We set an adaptation criterion with an error threshold of 0.02 on all
velocity components and $10^{-2}$ on the geometry. */

event adapt (i++) {
  astats s = adapt_wavelet ({cs,u}, (double[]){1e-2,0.02,0.02,0.02}, // fixme: these are not seen....
			    maxlevel, 4);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
