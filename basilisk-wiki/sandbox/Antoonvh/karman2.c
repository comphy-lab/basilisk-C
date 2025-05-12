/**
# Bénard–von Kármán Vortex Street for flow around a cylinder with particle tracers

[*The* example](/src/examples/karman.c) of 2D viscous flow around a
simple solid boundary: Fluid is injected to the left of a channel
bounded by solid walls with a slip boundary condition. Passive tracer
particles are injected near the cylinder.

![Animation of the tracer field](karman2/f.mp4)(loop)

We use the centered Navier-Stokes solver, with embedded boundaries and
advect the passive tracer particles `f`. and `f2` */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

Particles f, f2;

double Reynolds = 150.;
int maxlevel = 9;
face vector muv[];

/**
The domain is 5 units long, centered vertically. */

int main() {
  L0 = 5.;
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;
  run(); 
}

/**
We set a constant viscosity corresponding to a Reynolds number of 150,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*0.125/Reynolds;
  boundary ((scalar*){muv});
}
/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);

event init (t = 0) {
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
  f  = new_tracer_particles(0);
  f2 = new_tracer_particles(0);
 
  /**
  We set the initial velocity field. */
  foreach()
    u.x[] = cs[] ? 1. : 0.;
}
/**
Particles are added to the lists...
 */
event add_particles (t += 0.1) {
  double R = 0.07, A = 0.4;
  for (double thet = 0; thet < 2*pi; thet += pi/41) {
    particle p = {.x = R*sin(thet), .y = R*cos(thet)};
    set_a_particle_attributes (&p);
    add_particle (p, f);

    p.x -= A;
    set_a_particle_attributes (&p);
    add_particle (p, f2);
  }
}
/**
Outflow for the tracer particles is implemented with the following two
functions:
 */
void wrap_remove_particles (Particles p) {
    remove_particles (p , x > 4.4);
}

event outflow (i++) {
  foreach_P_in_list(tracer_particles) {
    wrap_remove_particles (P);
  }
}
/**
We adapt according to the error on the embedded geometry and velocity
fields. */
event adapt (i++) {
  adapt_wavelet ({cs, u}, (double[]){1e-3, 3e-2, 3e-2}, maxlevel, 4);
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We produce animations of the tracer field... */

event movies (i += 3; t <= 15) {
  view (fov = 5.5, width = 1024, height = 280, tx = -0.4);
  translate (z = 0.05)
    draw_vof("cs", "fs", filled = -1, fc = {0.4, 0.4, 0.4});
  scatter (f);
  scatter (f2, pc = {0.4, 0.6, 0.2}, s = 10);
  save ("f.mp4");
}


