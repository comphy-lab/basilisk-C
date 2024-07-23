/**
# Decaying two-dimensional turbulence

We solve the two-dimensional incompressible Euler equations using a
vorticity--streamfunction formulation. */

#include "grid/multigrid.h"
#include "navier-stokes/stream.h"

/**
The domain is square of size unity by default. The resolution is
constant at $256^2$. */

int main() {
  init_grid (256);
  run();
}

/**
The initial condition for vorticity is just a white noise in the range
$[-1:1]$ .*/

event init (i = 0) {
  double a = 1. [0,-1];
  foreach()
    omega[] = a*noise();
}

/**
We generate images of the vorticity field every 4 timesteps up to
$t=1000$. We fix the colorscale to $[-0.3:0.3]$.

![Evolution of the vorticity](turbulence/omega.mp4)(autoplay loop) */

event output (i += 4; t <= 1000) {
  output_ppm (omega, min = -0.3, max = 0.3, file = "omega.mp4");
}
