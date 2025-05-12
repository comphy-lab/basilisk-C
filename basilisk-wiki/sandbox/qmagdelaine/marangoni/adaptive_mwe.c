/**
# Example for the integral formulation of the surface tension 

We show here an example using an integral formulation of the surface tension
force, as depicted in
[Abu-Al-Saud, Popinet, Tchelepi, 2018](https://hal.archives-ouvertes.fr/hal-01706565/),
adapted to VOF.

For the moment, a lot of limitations remain. Nevertheless, this example already
shows that the current implementation works on an adaptive mesh (or at least seem
to work, we don't compare and quantitative data here).

![Pressure field](adaptive_mwe/pressure.mp4)

![Horizontal velocity field](adaptive_mwe/u_x.mp4)

*/

#define VIDEO 1

/**
## Includes */

#include "fractions.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"

#include "integral_vof_quasi1D.h"

#include "view.h"
#define BG 0.8 // light gray for background
#define DG 0. // dark gray

/**
## Paremeters */

#define L 4. // size of the box
#define LEVEL 7
#define MIN_LEVEL 5
#define T_END 3.
#define DELTA_T (T_END/100.)

#define h0 1.
#define dh 0.35 // 0.1
#define liquid_theta (45.*pi/180.)

/**
We define the shape of the interface. */

#define shape(x, y) (h0 + dh*cos(3.*pi*x/L) - y)

/**
We define a function to set and update the values of the surface tension on the
interface. We also define a scalar field to store it. */

void surface_tension (scalar gamma) {
  #if TREE
    gamma.refine = gamma.prolongation = curvature_prolongation;
    gamma.restriction = curvature_restriction;
  #endif
  foreach()
    gamma[] = (interfacial(point, f) ? 1. : nodata);
  boundary({gamma});
}
scalar liquid_sigma[];

/**
As with [contact.h](/src/contact.h), we define vector field to store the height
functions assciated to the VOF tracer f. */

vector h[];

int main() {
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);
  
  NITERMIN = 2;

  mu1 = 1.;
  mu2 = mu1/1.;

  /**
  We associate the height functions $h$ to the VOF tracer $f$. */

  f.height = h;
  run();
}

/**
## Boundary conditions

We set a free flow boundary conditions on top. */

u.n[top] = neumann(0.);
uf.n[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);

/**
## Initialization

In the init event, we initialize the surface tension. */

event init (i = 0) {
  fraction (f, shape(x, y));
  boundary({f});
  surface_tension (liquid_sigma);
  f.sigma = liquid_sigma;
}

/**
## Surface tension update

Since the surface tension is defined only on the interface (in view of having 
a non-uniform surface tension), we have to update it at each time step. */

event acceleration (i++) {
  surface_tension (liquid_sigma);
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){1e-4, 1e-1, 1e-1}, LEVEL, MIN_LEVEL);
}
#endif

/**
## Movie and outputs */

#if VIDEO
event movie_view (t = 0.; t += DELTA_T; t <= T_END) {
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  char legend[100];
  sprintf(legend, "t = %0.2g", t);
  int ratio = 2;
  
  view (fov = 19.1/ratio, width = ratio*320, height = 320, samples = 1,
        bg = {BG, BG, BG}, tx = - (0.5 + X0/L0),
        ty = - (0.5/ratio + Y0/L0), relative = false);
  clear();
  cells();
  draw_vof("f", edges = true, lw = 2., lc = {DG, DG, DG}, filled = 0);
  squares ("p", min = -2., max = 2., linear = false, map = cool_warm);
  draw_string(legend, 1, size = 30., lw = 2.);
  save ("pressure.mp4");
  
  view (fov = 19.1/ratio, width = ratio*320, height = 320, samples = 1,
        bg = {BG, BG, BG}, tx = - (0.5 + X0/L0),
        ty = - (0.5/ratio + Y0/L0), relative = false);
  clear();
  cells();
  draw_vof("f", edges = true, lw = 2., lc = {DG, DG, DG}, filled = 0);
  squares ("u.x", min = -0.1, max = 0.1, linear = false, map = cool_warm);
  draw_string(legend, 1, size = 30., lw = 2.);
  save ("u_x.mp4");
}
#endif
