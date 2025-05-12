/**
# Unsteady flow visualization based on the Bénard–von Kármán Vortex Street

[*The* example](/src/examples/karman.c) of 2D flow around a cylinder, Fluid is injected at the left of a channel formed by
free slip walls. Passive tracer particles
are injected near the cylinder.

We visualize the concepts of path lines, streak lines and stream lines:

![Path line (blue) and streak line (green). For the second path line,
 it is clear to see how a path line is the history of a fluid
 element](karman3/path_streak.mp4)

Furthermore, we  show the  streamfunction (in  the channel's  frame of
reference), vorticity field and a bunch of tracers.

![Stream lines, particles and vorticity field (color). Indeed, the
 unsteady streamlines do not represent material lines and cross particle
 paths](karman3/abunch_omg_psi.mp4)
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

Particles abunch, streak, path;

double Reynolds = 150.;
int maxlevel = 9;
face vector muv[];

/**
The domain is 6 units long, centered vertically. */


long int pcounter = 0, tagcounter = 0;

int main() {
  L0 = 6.;
  origin (-.5, -L0/2.);
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

u.n[left]  = dirichlet(cs[] ? 1.: 0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = dirichlet(cs[] ? 1.: 0);
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
  /**
     A bunch of particles are randomly added in the domain.
  */
  int np = 200, iadded = 0;
  abunch  = new_tracer_particles(0);
  while (iadded < np) {
    double xp = X0 + fabs(noise())*L0;
    double yp = Y0 + fabs(noise())*L0;
    if (interpolate (cs, xp, yp) > 0.99) {
      particle pnew = {.x = xp, .y = yp};
      set_a_particle_attributes (&pnew);
      add_particle (pnew, abunch);
      iadded++;
    }
  }
  /**
The streak and path line start out as empty lists.
   */
  streak = new_tracer_particles(0);
  path = new_particles(0);
  /**
  We set the initial velocity field. */
  foreach()
    u.x[] = cs[] ? 1. : 0.;
}
/**
Particles are added to the lists...
 */
event add_particles (t += 0.005) {
  particle p = {.x = X0, .y = 0.05, .tag = pcounter++};
  set_a_particle_attributes (&p);
  add_particle (p, streak);

  foreach_particle_in (streak) {
    if (p().tag == tagcounter) {
      particle pnew = {.x = x, .y = y};
      set_a_particle_attributes (&pnew);
      add_particle (pnew, path);
      break;
    }
  }
}
/**
## Outflow of particles
 */

event outflow (i++) {
  remove_particles (streak , x > 5.41);
  int pnb = pn[path]; //Check if the path-line head is removed.
  remove_particles (path , x > 5.4);
  if (pn[path] < pnb) { 
    remove_particles (path, 1); // Remove all pathline particles
    tagcounter = pcounter; // Start following a new one from the streakline list.
  }
}
/**
We adapt according to the estimated error on the embedded geometry and
velocity fields. */
event adapt (i++) {
  adapt_wavelet ({cs, u}, (double[]){1e-4, 1e-2, 1e-2}, maxlevel, 7);
}

/**
We log the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
## Movie maker event

   We produce animations of our particles, vorticity and stream function.
*/
scalar psi[]; // Remember stream function as starting guess for iterative solver
event movies (t += .025; t <= 17) {
  view (fov = 5., width = 1000, height = 250, tx = -0.42, ty = -0.02);

  /**
     The stream function ($\psi$) is computed as a function of the
     vorticity field ($\omega$) following their (Poisson) relation.

     $$\nabla^2 \psi = -\omega$$
     
   */
  scalar omg[];
  foreach()
    omg[] = 0;
  vertex scalar psif[];
  if (i > 5) {
    vorticity (u, omg);
    foreach() {
      psi[] = 0;
      if (cs[] > 0 && cs[] < 1) {
	coord n, p;
	embed_geometry (point, &p, &n);
	omg[] = embed_vorticity (point, u, p, n);
      }
      if (cs[] == 0)
	omg[] = 0;
    }
    /**
The boundary conditions for $\psi$ are determined by our frame of
reference. Furthermore, note that the condition on the cylinder is
only approximate.
     */
    psi[embed] = fabs(y) > 0.25 ? dirichlet(y < 0 ? 0.5 : -0.5) : dirichlet(0.);
    poisson (psi, omg, alpha = fs, tolerance = 1e-2);
    // Interpolate to vertices:
    foreach() {
      if (y > 0.5)
	psi[] = -0.5;
      if (y < -0.5)
	psi[] = 0.5;
    }
    foreach_vertex()
      psif[] = (psi[0,0] + psi[-1,0] + psi[-1,-1] + psi[0,-1])/4.;
    

    // Combined tracers+omg+psi  
    translate (z = 0.02)
      for (double val = -.3; val <= .301; val += 0.075) 
	isoline ("psif", val, lw = 2);
  }
  translate (z = 0.05)
    draw_vof("cs", "fs", filled = -1, fc = {0.6, 0.6, 0.6});
  scatter (abunch);
  squares ("omg", min = -5, max = 5, map = blue_white_red, linear = true);
  draw_string (" Tracers, streamlines and vorticity field", size = 43,  pos = 1, lw = 2);
  save ("abunch_omg_psi.mp4");
  save ("abunch.png");
  // Combined Steak+path
  clear();
  translate (z = 0.05)
    draw_vof("cs", "fs", filled = -1, fc = {0.6, 0.6, 0.6});
  scatter (streak, pc = {0.4, 0.6, 0.2}, s = 4);
  scatter (path, pc = {0.2, 0.4, 0.6}, s = 3);
  draw_string (" Pathline (blue)", lc = {0.1, 0.2, 0.3}, lw = 2, pos = 1);
  draw_string ("streak line (green) ", lc = {0.2, 0.3, 0.1}, lw = 2, pos = 2);
  save ("path_streak.mp4");

}


