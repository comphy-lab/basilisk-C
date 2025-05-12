/**
# A dipolar vortex in a large forest

![Mind the dynamically reconstructed geometry.](test_reconst/o.mp4)
 */
#include "embed.h"

double yp = 0; 
double forest (double x, double y) {
  if (y > yp ) 
    return 1.;
  else
    return sin(pi*x) + cos(pi*y) + 0.6;
}

#define GEOM(x,y) forest (x, y)

#include "reconstruct.h"
#include "navier-stokes/centered.h"
#include "view.h"

face vector muc[];
double Re = 500.;
u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

int main() {
  L0 = 35;
  X0 = Y0 = -L0/2.;
  mu = muc;
  run();
}

double xo = 0.1, yo = 10;
#define RAD (sqrt(sq(x - xo) + sq(y - yo)))
#define ST (-(x - xo)/RAD)
event init (i = 0) {
  scalar psi[];
  double k = 3.83170597;
  refine (level <= 6); //Init GEOM
  refine (RAD < 2.0 && level <= 9);
  refine (RAD < 1.0 && level <= 10);
  foreach() 
    psi[] = ((RAD > 1)*(ST/RAD) +
	     (RAD < 1)*((-2*j1(k*RAD)*ST/(k*j0(k))) + (RAD*ST)));
  boundary ({psi});
  foreach() {
    u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta));
    u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta);
  }
  boundary (all);
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]/Re;
  boundary ((scalar*){muc});
}

event movie (t += 0.1; t <= 25) {
  scalar omega[];
  vorticity (u,omega);
  view (fov = 21, width = 950, height = 500);
  translate (x = -L0/1.95) {
    box (notics = true);
    draw_vof ("cs", "fs");
    draw_vof ("cs", "fs", filled = -1, fc = {0.8,0.1,0.8});
    squares ("omega", map = cool_warm, min = -5, max = 5);
  }
  translate (x = L0/1.95)
    cells();
  save ("o.mp4");
}

event adapt (i++) {
  adapt_wavelet ({u.x, u.y}, (double[]){0.01, 0.01}, 10);
}
  
