/**
#
On Ubuntu 20.04 LTS, when I run file examples/karman.c,
I don't see any problem with the file f.mp4, but
the file vort.mp4 looks differently from the video on this website.

Before reporting this, I updated both Ubuntu and Basilisk.

Near the top border I see (dark) red and
near the bottom border I see (dark) blue.
The width of these regions is comparable to the cylinder diameter.
In the middle of the geometry, as the flow advances,
there might not be any change.

# Bénard–von Kármán Vortex Street for flow around a cylinder at Re=160

![Animation of the vorticity field.](bug_karman/vort.mp4)(loop)
 
*/

#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
double Reynolds = 160.;
int maxlevel = 9;
face vector muv[];

/**
The domain is eight units long, centered vertically. */

int main() {
  L0 = 8.;
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;

  /**
  When using bview we can interactively control the Reynolds number
  and maximum level of refinement. */
  
  display_control (Reynolds, 10, 1000);
  display_control (maxlevel, 6, 12);
  
  run(); 
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*0.125/Reynolds;
}

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);

event init (t = 0)
{

  /**
  The domain is the intersection of a channel of width unity and a
  circle of diameter 0.125. */

  vertex scalar phi[];
  foreach_vertex() {
    phi[] = intersection (0.5 - y, 0.5 + y);
    phi[] = intersection (phi[], sq(x) + sq(y) - sq(0.125/2.));
  }
  boundary ({phi});
  fractions (phi, cs, fs);

  /**
  We set the initial velocity field. */
  
  foreach()
    u.x[] = cs[] ? 1. : 0.;
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We produce animations of the vorticity and tracer fields... */

event movies (i += 4; t <= 15.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  boundary ({m});
  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true, mask = m);
  output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      linear = false, min = 0, max = 1, mask = m);
}

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
}
