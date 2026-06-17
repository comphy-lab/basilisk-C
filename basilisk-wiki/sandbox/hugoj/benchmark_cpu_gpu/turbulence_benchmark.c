/**
# Decaying two-dimensional turbulence

This code is a modified version of [turbulence.c](https://basilisk.fr/src/examples/turbulence.c).
I added output at beginning and end of simulation, and when setting BENCHMARK=1
the movie is not generated and the GPU is not outputing anything to the screen. 

TODO:
- fixme: the movie generated stutters when runnning with the movie output on HPC
GPU. Temporary fix: outputs only snapshots
*/

#include "navier-stokes/stream.h"

/**
The domain is square of size unity by default. The default resolution
is constant at $256^2$ but can be change using command line
arguments. */

int main (int argc, char * argv[]) {
  init_grid (argc > 1 ? atoi(argv[1]) : 256);
  run();
}

/**
The initial condition for vorticity is just a white noise in the range
$[-1:1]$ .*/

event init (i = 0) {
  double a = 1. [0,-1];
  foreach (cpu)
    omega[] = a*noise();
}

/**
We generate images of the vorticity field every 4 timesteps up to
$t=1000$. We fix the colorscale to $[-0.3:0.3]$.

The expected output is like the movie below

![Evolution of the vorticity](https://basilisk.fr/src/examples/turbulence/omega.mp4)(autoplay loop) */

#if !BENCHMARK
event output (t += 1; t <= 1000) {
  output_ppm (omega, min = -0.3, max = 0.3, file = "omega.mp4");
}
#endif

event snap_ini (t=0) {
  output_ppm (omega, min = -0.3, max = 0.3, file = "omega_ini");
}
event snap_end (t=1000) {
  output_ppm (omega, min = -0.3, max = 0.3, file = "omega_end");
}



#if !BENCHMARK
#if _GPU && SHOW
event display (i++)
  output_ppm (omega, min = -0.3, max = 0.3, fps = 30, fp = NULL);
#endif
#endif
/**

[Benchmark on GPUs](/src/grid/gpu/Benchmarks.md#two-dimensional-turbulence)
*/
