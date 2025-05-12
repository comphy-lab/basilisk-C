/**
# Falling Droplet in Moving Reference Frame

We want to simulate a falling droplet in a reference frame
centered on the center of the droplet. In this situation,
the droplet falls without leaving the domain and remaining
fixed to the position of the initialization.
The implementation is very naive, and it consists in a
feedback control of the inlet velocity, which is used to
compensate the gravitational acceleration. Better solutions
and suggestions are welcome.

![Evolution of the droplet position and the vorticity field](fallingdropmrf/movie.mp4)(width="600" height="600")

The droplet starts to fall from a zero initial velocity and
to experience the sudden increase in the inlet velocity. Therefore,
it takes some time to stabilize. Finally, it reaches a quite
steady behavior, although it starts to slighly shift along the
x axis.
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"

/**
We initialize the additionals fields required for this simulation:
*fx* and *fy* are used to compute the centroid of the droplet;
*omega* is the vorticity field. */

scalar fx[], fy[];
scalar tr[];

/**
We use inlet boundary condition on the bottom of the domain,
with a velocity that varies depending on the acceleration of
the droplet, and outflow boundary conditions on the top. */

double uin;
u.n[bottom] = dirichlet(uin);
u.t[bottom] = dirichlet(0);
p[bottom] = neumann(0);

u.n[top] = neumann(0);
u.t[top] = neumann(0);
p[top] = dirichlet(0);

/**
We set the level of refinement, the total simulation
time, the initial radius of the droplet, and we declare
the variables to store the centroid of the droplet, and
the velocity of that centroid. */

int maxlevel, minlevel;
double tEnd = 0.05;
double R0 = 200.e-6;
double xcentroid, ycentroid;
double xcentroid_old, ycentroid_old;
double uxcentroid, uxcentroid_old;
double uycentroid, uycentroid_old;

int main (void)
{
  /**
  We set the physical properties for the simulation. */

  rho1 = 800.; rho2 = 5.;
  mu1 = 1.138e-3; mu2 = 1.78e-5;

  f.sigma = 0.03;

  /**
  We define the problem geometry and we include the
  gravitational acceleration using the [reduced.h](/src/reduced.h)
  module. */

  minlevel = 2;
  maxlevel = 8;

  L0 = 6.*(6.*R0);
  G.y = -9.81;

  size (L0);
  origin (0., 0.);
  init_grid (1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x - 0.5*L0) - sq(y - 0.3*L0))

event init (i = 0) {
  /**
  We initialize the droplet and we compute the position
  of its centroid. */

  fraction (f, circle (x, y, R0));
  foreach() {
    fx[] = f[]*x;
    fy[] = f[]*y;
  }
  xcentroid = statsf(fx).sum / statsf(f).sum;
  ycentroid = statsf(fy).sum / statsf(f).sum;
}

/**
In this event, the centroid position is computed,
as well as the inlet velocity. The shifting
of the position of the droplet centroid with
respect to the target (initial) position is used
to compute a velocity which is imposed with opposite
sign at the bottom of the domain. */

event bcs (i++) {
  foreach() {
    fx[] = f[]*x;
    fy[] = f[]*y;
  }
  xcentroid_old = xcentroid;
  ycentroid_old = ycentroid;
  xcentroid = statsf(fx).sum / statsf(f).sum;
  ycentroid = statsf(fy).sum / statsf(f).sum;

  uin = -(ycentroid - 0.3*L0)/dt;
  u.n[bottom] = dirichlet (uin);
}

/**
We refine the mesh according to the volume fraction and the
velocity field. */

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1.e-3,1.e-3,1.e-3,1.e-3}, maxlevel, minlevel);
}

/**
The following events are for post-processing purposes:
print the inlet velocity value, compute the vorticity field,
write a video with the evolution of the interface and the
vorticity field. */

event outputfile (i++) {
  fprintf (stdout, "%f %f\n", t, uin);
}

event movie (t += 0.0001; t <= tEnd) {
  scalar omega[];
  vorticity (u, omega);

  clear();
  box();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("omega", linear = false,
           min = -6000., max = 6000.,
           map = blue_white_red);
  save ("movie.mp4");
}
