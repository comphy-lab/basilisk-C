/**
# A source of waves.
   
   ![Waves originating from a central source in a square container](source/p.mp4)

   ![Time-averaged intensity reveals the interference pattern](source/P.mp4)
*/

#include "wave.h"

scalar P[];

int main() {
  L0 = 100;
  X0 = Y0 = -L0/2.;
  N = 512;
  DT = 0.1;
  run();
}

event acceleration (i++) {
  foreach() 
    a[] = sin(t)*exp(-sq(x) - sq(y));
}

event movies (t += 0.1) {
  output_ppm (s, file = "p.mp4", n = 300, min = -1, max = 1);
}

event intensity (t += 0.1) {
  foreach()
    P[] += sq(s[]);
  output_ppm (P, file = "P.mp4", n = 300);
}

event stop (t = 150);
