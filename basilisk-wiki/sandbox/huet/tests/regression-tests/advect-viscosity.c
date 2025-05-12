/**
# Advection of a spherical capsule featuring a viscosity jump

In this file we test the advection of a capsule featuring a non-unity
viscosity ratio, across a perodic boundary. This file tests the 
functions in [viscosity-ft.h](../../src/lagrangian_caps/viscosity-ft.h) as well
as the lag2eul function in [ibm-ft/h](../../src/lagrangian_caps/ibm-ft.h). */

#define LEVEL 6
#define LAG_LEVEL 4
#define RADIUS .2
#define L0 1.
#define T_END 1.

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/viscosity-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"

int main(int argc, char* argv[]) {
  origin(-.5*L0, -.5*L0, -.5*L0);
  N = 1 << LEVEL;
  init_grid(N);
  periodic(left);
  TOLERANCE = HUGE;
  stokes = true;
  DT = 1.e-2;
  run();
}

event init (i = 0) {
  activate_spherical_capsule(&CAPS(0), level = LAG_LEVEL, radius = RADIUS);
}

event impose_u (i++) {
  foreach() {
    double x0 = x + L0 -t*L0*T_END;
    double y0 = y + L0;
    double a = -2*sq(sin(pi*x0))*sin(pi*y0)*cos(pi*y0)*cos(pi*t/T_END);
    double b = -2*sin(pi*x0)*cos(pi*x0)*sq(cos(pi*y0))*cos(pi*t/T_END);
    double theta = pi/4.;
    u.x[] = a*cos(theta)/2 + L0*T_END;
    u.y[] = b/2;
    u.z[] = a*sin(theta)/2;
  }
}

event adapt (i++) {
  tag_ibm_stencils(&CAPS(0));
  adapt_wavelet({stencils}, (double []){1.e-2}, maxlevel = LEVEL);
  generate_lag_stencils(&CAPS(0));
}

event progress (i++) {
  fprintf(stdout, "i=%d, t=%g\n", i, t);
  fflush(stdout);
}

/** We compute the time evolutions of the normalized area and volume 
of the capsule */
event movie (i++) {
  view(fov = 18.9, bg = {1,1,1}, theta = 5*pi/6, psi = 0., phi = pi/8);
  clear();
  draw_lag(&CAPS(0), lw = .5, edges = true, facets = false);
  cells(n = {0,0,1});
  squares("I", n = {0,0,1}, min=-.5, max=1.5);
  save("viscosity.mp4");
}

/** At the end of the simulation, we also compare the position of the
membrane to its initial position. With an asymptotically fine space
and time resolutions, they would coincide. */
event output (t = T_END/2) {
  foreach() 
    if (I[] > 1.e-10) 
      fprintf(stderr, "%.4g, %.4g, %.4g, %.4g\n", x, y, z, I[]);
}

event end (t = T_END/2) {
  return 0;
}

/**
## Results

![Viscosity field advected through a periodic boundary](advect-viscosity/viscosity.mp4)
*/