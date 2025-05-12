/**
# Drop impact on a curved surface.
I made this simulation for testing the vof+embed+reduced system.

![Vorticity field and gas-liquid interface](embedfalldrop/movie.mp4)(width="800" height="600")
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "tension.h"
#include "view.h"

/**
We initialize the vorticity field and the maximum and minimum
levels of refinement. */

int maxlevel = 8;
int minlevel = 5;

/**
We define the radius of the solid obstable, and the total
simulation time. */

double R0 = 0.4e-2;
double tEnd = 1.e-1;

/**
We set no-slip conditions on the solid obstacle (In theory
more realistic boundary conditions than symmetry everywhere
should be set for the other boundaries). */

u.n[embed] = dirichlet (0.); 
u.t[embed] = dirichlet (0.);
p[embed] = neumann (0.);

int main (void) {

  /**
  We initialize the physical properties of the
  two-phase system and the gravity value. */

  rho1 = 800.; rho2 = 5.;
  mu1 = 1.e-3; mu2 = 1.e-5;

  f.sigma = 0.073;
  G.y = -9.81;

  L0 = 1.e-2;
  origin (0., 0.);
  init_grid (1 << maxlevel);
  run();
}

event init (i = 0) {
  
  /**
  We initialize the field *cs* which defines the region where
  the fluid dynamics is resolved: $cs > 0$. */

  solid (cs, fs, -(sq(R0) - sq(x) - sq(y)));

  /**
  We initialize a liquid droplet at coordinates
  $(x_c,y_c) = (0.2L_0, 0.7L_0)$ and radius $R = 0.1L_0$. */

  fraction (f, -(sq(x - 0.2*L0) + sq(y - 0.7*L0) - sq(0.1*L0)));
}

/**
We refine the region around the interface of the droplet. */

event adapt (i++) {
  adapt_wavelet ({cs,f,u}, (double []){1.e-3, 1.e-3, 1.e-3, 1.e-3}, maxlevel);
}

/**
The following events are for post-processing purposes:
compure the vorticity field, write a video with the
evolution of the interface and the vorticity field,
stop the simulation. */

event movie (t += 0.001; t <= tEnd) {
  scalar omega[];
  vorticity (u, omega);
  
  clear();
  box();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("cs", "fs", filled = -1, fc = {0.3,0.3,0.3});
  draw_vof ("f", lw = 1.5);
  squares ("omega", linear = false, 
           min = -1000., max = 1000., 
           map = blue_white_red);
  save ("movie.mp4");
}
