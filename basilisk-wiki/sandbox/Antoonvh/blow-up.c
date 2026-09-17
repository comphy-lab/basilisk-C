/**
# Blow-up(?) by vortex stretching

Celebrating the [news of the cracked millennium
problem](https://openai.com/index/navier-stokes-solution/), we model a
vortex that is forced to stretch.

![Azimuthal velocity field of a stretching vortex. It appears to be smiling](blow-up/w.mp4)

~~~gnuplot Maximum azimuthal velocity does not blow up
set xlabel 'time'
set ylabel 'Max. Azim. vel.'
set key off
set grid
plot 'out' w l lw 2 t 'data'
~~~
 */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/swirl.h"
#define MUV (5e-4)
face vector av[];

const face vector muc[] = {MUV, MUV};

int min_level = 7, max_level = 11;

int main() {
  mu = muc;
  a = av;
  L0 = 12;
  X0 = -L0/2;
  N = 256;
  run();
}

event init (t = 0) {
  double vk = 1;
  foreach()
    w[] = vk*y*exp(-sq(x) - sq(y));
}
/**
## Forcing
   
We have a stretching force near the vortex centre ...
*/
event acceleration (i++) {
  double ak = 1;
  foreach_face(x)
    av.x[] = x*ak*exp(-sq(x) - sq(y));
}
/**
  ... and a swirling force.
*/
event tracer_diffusion(i++) {
  double ak = 0.1;
  foreach()
    w[] += ak*dt*y*exp(-sq(x) - sq(y));
}

event movie (t += 0.1) {
  output_ppm (w, file = "w.mp4", n = 300, min = 0, max = 1);
  output_ppm (u.x, file = "ux.mp4", n = 300, min = -1, max = 1);
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (lev, file = "lev.mp4", n = 300, min = min_level - 0.5, max = max_level + 0.5);
}

event statistics_event (t += 0.05) {
  printf ("%g %g %ld\n", t, statsf(w).max, grid->tn);
}

event adapt (i++) {
  adapt_wavelet ({u, w},(double[]){1e-3, 1e-3, 1e-3}, max_level, min_level);
}

event stop (t = 10) {
  return 1;
}
