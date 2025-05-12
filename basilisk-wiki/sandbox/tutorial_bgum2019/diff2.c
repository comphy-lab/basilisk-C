/**
# Tutorial: Using the new diffusion solver

![Results from "our" solver](diff2/s.mp4)

![Its is compatible with adaptivity too!](diff2/level.mp4)

We closely follow the [previous](diff.c) [examples](diff_v.c) from
this tutorial. 
*/

#include "diffusion2.h"
#include "run.h"

scalar s[];

/**
The value of the time-stepping parameter `DT`, reveals the main
advantage of the [diffusion](/src/diffusion.h) solver that comes with
Basilisk. Our forward-in-time integrators are prone to numerical
instabilities and require a small time step for the succes of the
simulation.
*/
int main() {
  DT = 0.005; //ten times smaller than the implicit solver
  run();
}

event init (t = 0) {
  foreach() 
    s[] = exp(-(sq(x - 0.5) + sq(y - 0.5))*10.); 
}

event mov (t += 0.1) {
  output_ppm (s, file = "s.mp4", n = 256, min = -1, max = 1);
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (lev, file = "level.mp4", n = 256, max = 6);
}

event diff (i++) {
  dt = dtnext (DT);
  const face vector kap[] = {0.01, 0.01};
  diffusion_midpoint (s, dt, kap);
}

event lot (i += 5) {
  static FILE * fp = fopen ("data", "w");
  fprintf (fp, "%g %d %.8g\n", t, i, statsf(s).sum);
}

event adapt (i++)
  adapt_wavelet ({s}, (double[]){0.01}, 6);

event stop (t = 10);

