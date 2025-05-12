/**
#Spurious currents on a naively adapted grid.

This case concerns a bubble with turface tension

![animation of the horizontal velocity. Adaptivity is switched off](spurcur/ux2.mp4)

 */
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

int main(){
  f.sigma = 1;
  TOLERANCE = 1e-6;
  init_grid (1 << 8);
  run();
}

event init (t = 0) {
    fraction (f, sq(x -0.5 ) + sq(y - 0.5) - sq(0.25)); 
}

event adapt(i = 500; i++)
  adapt_wavelet({f}, (double[]){0.001}, 9);

event movie (i += 5)
  output_ppm(u.x, file = "ux2.mp4", n = 500,
             min = -0.02, max = 0.02);

event stop (i = 1000);