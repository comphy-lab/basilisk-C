/**
# A Wave-equation solver

![result for a bump](bump_wave/s.mp4)
 */
#include "grid/multigrid.h"
#include "wave.h"

int main() {
  L0 = 20;
  X0 = Y0 = -L0/2.;
  N = 256;
  periodic (left);
  periodic (bottom);
  DT = 0.5;
  TOLERANCE = 1e-7;
  run();
}

event init (t = 0) {
  foreach() {
    s[] = exp(-sq(x)-sq(y));
  }
}

event outputer (i++) {
  output_ppm (s, file = "s.mp4", n = 300, min = -1, max = 1);
}

event stop (t = 100);
