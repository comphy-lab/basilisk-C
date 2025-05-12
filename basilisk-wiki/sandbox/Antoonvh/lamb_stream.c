/**
# A test for the $\omega-\psi$ Navier-Stokes solver

We advect an initial Lamb-Chaplygin vortex dipole, subject to a viscous force.

![This looks OK](lamb_stream/movie.mp4)

The result may be compared against he one [obtained with the `mac.h`](lamb_mac.c)
or the [`centered.h`](lamb_centered.c) solver.
*/

#include "grid/multigrid.h"
#include "navier-stokes/stream.h"
#include "diffusion.h"

#define RAD (sqrt(sq(x) + sq(y)))
#define ST (x/RAD)

const face vector nu[] = {0.01, 0.01};

int main() {
  L0 = 20;
  X0 = Y0 = -L0/2.;
  N = 256;
  run();
}

event init (t = 0) {
  double k = 3.83170597;
  foreach()
    omega[] = -sq(k)*(RAD < 1)*(-2*j1(k*RAD)*ST/(k*j0(k)));
}
/**
Despite the equation listed under `stream.h`, it does not include viscous diffusion of vorticity.
*/
event tracer_diffusion (i++)
  diffusion (omega, dt, nu);
  
event movie (t += 0.1; t < 10)
  output_ppm (omega, file = "movie.mp4", n = 400);