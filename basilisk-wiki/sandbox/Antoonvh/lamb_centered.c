/**
# A test for the centered Navier-Stokes solver

We advect an initial Lamb-Chaplygin vortex dipole, subject to a viscous force.

![This looks OK](lamb_centered/movie.mp4)

The result may be compared against he one [obtained with tme `mac.h`](lamb_mac.c) solver.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

#define RAD (sqrt(sq(x) + sq(y)))
#define ST (x/RAD)

const face vector nu[] = {0.01, 0.01};

int main() {
  mu = nu;
  L0 = 20;
  X0 = Y0 = -L0/2.;
  N = 256;
  run();
}

event init (t = 0) {
  double k = 3.83170597;
  scalar psi[];
  foreach()
    psi[] = ((RAD > 1)*((1/RAD))*ST +
	     (RAD < 1)*((-2*j1(k*RAD)*ST/(k*j0(k))) + (RAD*ST)));
  boundary({psi});
  foreach() {
    u.x[] = -(psi[0, 1] - psi[0, -1])/(2*Delta);
    u.y[] = (psi[1] - psi[-1])/(2*Delta);
  }
}

event movie (t += 0.1; t < 10){
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, file = "movie.mp4", n = 400);
}
