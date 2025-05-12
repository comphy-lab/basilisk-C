/**
# Example for the integral formulation: Marangoni and contact angles

We show here an example using an integral formulation of the surface tension
force, as depicted in
[Abu-Al-Saud, Popinet, Tchelepi, 2018](https://hal.archives-ouvertes.fr/hal-01706565/),
adapted to VOF.

For the moment, a lot of limitations remain. Nevertheless, this example already
shows that this new implementation of the surface tension allow to describe
contact angles and Marangoni flows.

![Velocity in the $x$ direction](marangoni_test/u_x.mp4)

The velocty field seems nice, we see the Marangoni flows induced by the surface
tension gradient and the respect of the wetting condition.

![Pressure](marangoni_test/pressure.mp4)

In the pressure, we see some glitches at the contact line (on the right and on
the left) and some global fluctuations. These glitches appear when I try to set
the contact angle. */

#define VIDEO 1

/**
## Includes */

#include "fractions.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"

/**
Instead of using [tension.h](/src/tension.h) and [contact.h](contact.h), to
account for the surface tension and contact angles, we use
[integral_vof_quasi1D.h](integral_vof_quasi1D.h). */

#include "integral_vof_quasi1D.h"

#include "view.h"
#define BG 0.8 // light gray for background
#define DG 0. // dark gray

/**
## Paremeters */

#define L 4. // size of the box
#define LEVEL 7
#define T_END 10.
#define DELTA_T (T_END/100.)
#define F_ERR 1e-10

/**
We set the value of the contact angle of the liquid on the walls. */

#define liquid_theta (45.*pi/180.)

#define h0 1.
#define dh 0.2
#define x0 (-asin(dh*6.*pi/L)*L/(6.*pi))
#define shape(x, y) (h0 + dh*cos(6.*pi*(x - x0)/L) - y)

/**
We define a function to set and update the values of the surface tension on the
interface. We also define a scalar field to store it. */

void surface_tension (scalar gamma) {
  #if TREE
    gamma.refine = gamma.prolongation = curvature_prolongation;
    gamma.restriction = curvature_restriction;
  #endif
  foreach()
    gamma[] = (interfacial(point, f) ? 1. - 0.5*cos(pi*(x/L0 - 0.5)) : nodata);
  boundary({gamma});
}
scalar liquid_sigma[];

/**
As with [contact.h](/src/contact.h), we define vector field to store the height
functions assciated to the VOF tracer f. */

vector h[];
FILE * fp = NULL;

int main() {
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);

  mu1 = 1.;
  mu2 = mu1/10.;

  f.height = h;
  run();
}

/**
# Boundary conditions on the pressure and the velocity

We set a free flow condition on the top and solid conditions and the three
other bondaries. */

u.n[top] = neumann(0.);
uf.n[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);

u.t[bottom] = dirichlet(0.);
uf.t[bottom] = dirichlet(0.);
u.t[left] = dirichlet(0.);
uf.t[left] = dirichlet(0.);
u.t[right] = dirichlet(0.);
uf.t[right] = dirichlet(0.);

/**
## Initialization

In the init event, we initialize the surface tension and set different
attributes to the VOF tracer $f$:

* its potentially variable surface tension,
* the contact angle of the liquid on the different boundaries (for the moment,
it is only avaible on the left and on the right). */

event init (i = 0) {
  fraction (f, shape(x, y));
  boundary({f});

  surface_tension (liquid_sigma);
  f.sigma = liquid_sigma;

  f.theta_left = liquid_theta;
  f.theta_right = liquid_theta;
}

/**
## Surface tension update

Since the surface tension is defined only on the interface, we have to update it
at each time step to follow the interface. */

event acceleration (i++) {
  surface_tension (liquid_sigma);
}

/**
## Videos */

#if VIDEO
event movie (t = 0.; t += DELTA_T; t <= T_END) {

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  char legend[100];
  sprintf(legend, "t = %0.2g", t);
  
  /**
  We write a first video with the pressure field. */
  
  view (fov = 19., width = 640, height = 640, samples = 1,
        bg = {BG, BG, BG}, tx = - (0.5 + X0/L0),
        ty = - (0.5 + Y0/L0), relative = false);
  clear();
  draw_vof("f", edges = true, lw = 2., lc = {DG, DG, DG}, filled = 0);
  squares ("p", min = -2., max = 2., linear = false, map = cool_warm);
  draw_string(legend, 1, size = 30., lw = 2.);
  save ("pressure.mp4");
  
  /**
  We write a second video with the horizontal velocity field. */
  
  view (fov = 19., width = 640, height = 640, samples = 1,
        bg = {BG, BG, BG}, tx = - (0.5 + X0/L0),
        ty = - (0.5 + Y0/L0), relative = false);
  clear();
  draw_vof("f", edges = true, lw = 2., lc = {DG, DG, DG}, filled = 0);
  squares ("u.x", min = -0.1, max = 0.1, linear = false, map = cool_warm);
  draw_string(legend, 1, size = 30., lw = 2.);
  save ("u_x.mp4");
}
#endif
