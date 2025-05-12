/**
# A Wave-equation solver

![Result for a travelling sine wave](wave/s.mp4)
 */
#include "grid/multigrid1D.h"
#include "wave.h"

int main() {
  L0 = 10;
  N = 256;
  periodic (left);
  DT = 0.5;
  TOLERANCE = 1e-7;
  run();
}

event init (t = 0) {
  foreach() {
    s[] = sin(pi*x);
    ds[] = -cos(pi*x)*pi; 
  }
}

event outputer (i++) {
  output_ppm (s, file = "s.mp4", n = 300, min = -1, max = 1);
}

event stop (t = 100);
