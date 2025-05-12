/**
# A test for the Markers and Cells Navier-Stokes solver

We advect an initial Lamb-Chaplygin vortex dipole, subject to a viscous force.

![The results are not very convincing.](lamb_mac/movie.mp4)

Using more refined grids appears to reduce the issues.
*/

#include "grid/multigrid.h"
#include "navier-stokes/mac.h"

#define RAD (sqrt(sq(x) + sq(y)))
#define ST (x/RAD)

int main() {
  nu = 0.01;
  CFL = 0.5;
  TOLERANCE = 1e-5;
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
  /**
   This bit is taken from [`stream.h`](http://www.basilisk.fr/src/naver-stokes/stream.h)
   */
  coord f = {-1.,1.};
  foreach_face()
    u.x[] = f.x*(psi[0,1] + psi[-1,1] - psi[0,-1] - psi[-1,-1])/(4.*Delta);
}

event movie (t += 0.1; t < 10){
  scalar omega[];
  foreach()
    omega[] = (u.y[1] - u.y[] - u.x[0,1] + u.x[])/Delta;  
  output_ppm (omega, file = "movie.mp4", n = 400);
}
