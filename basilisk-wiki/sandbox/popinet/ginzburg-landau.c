/**
# The complex Ginzburg--Landau equation (on GPUs)

This is the GPU version of [this example](/src/examples/ginzburg-landau.c). With a resolution of 2048^2^ it runs fast on the RTX 4090:

![](ginzburg-landau/out.txt)

![Evolution of the imaginary part $A_i$](ginzburg-landau/Ai.mp4)(autoplay) */

#include "grid/gpu/multigrid.h"
#include "run.h"
#include "diffusion.h"

scalar Ar[], Ai[], A2[];

/**
In this example, we only consider the case when $\alpha=0$. */

double beta;

/**
The generic time loop needs a timestep. We will store the statistics
on the diffusion solvers in `mgd1` and `mgd2`. */

double dt;
mgstats mgd1, mgd2;

/**
## Parameters

We change the size of the domain `L0`. */

int main() {
  beta = 1.5;
  size (500);
  init_grid (2048);
  run();
  system ("mv out out.txt");
}

/**
## Initial conditions 

We use a white noise in $[-10^{-4}:10^{-4}]$ for both components. */

event init (i = 0) {
  foreach() {
    Ar[] = 1e-4*noise();
    Ai[] = 1e-4*noise();
  }
}

/**
## Time integration */

event integration (i++) {

  /**
  We first set the timestep according to the timing of upcoming
  events. We choose a maximum timestep of 0.05 which ensures the stability
  of the reactive terms for this example. */

  dt = dtnext (0.05);

  /**
  We compute $|A|^2$. */

  foreach()
    A2[] = sq(Ar[]) + sq(Ai[]);

  /**
  We use the diffusion solver (twice) to advance the system from $t$
  to $t+dt$. */

  scalar r[], lambda[];
  foreach() {
    r[] = A2[]*beta*Ai[];
    lambda[] = 1. - A2[];
  }
  mgd1 = diffusion (Ar, dt, r = r, beta = lambda);
  foreach() {
    r[] = - A2[]*beta*Ar[];
    lambda[] = 1. - A2[];
  }
  mgd2 = diffusion (Ai, dt, r = r, beta = lambda);
}

/**
## Outputs

Here we create MP4 animations for both components. The `spread`
parameter sets the color scale to $\pm$ twice the standard
deviation. */

event movies (t += 0.2; t <= 150) {
  fprintf (stderr, "%g %g\n", t, sqrt(normf(A2).max));

  output_ppm (Ai, spread = 2, linear = true, file = "Ai.mp4", n = 1024);
}
