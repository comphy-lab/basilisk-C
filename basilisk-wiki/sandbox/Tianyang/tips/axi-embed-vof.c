/**
# Update of cm & fm 

When trying to combine axi-embed-twoPhase modules, metric cm and fm should be updated in the initialization step.
*/

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"

#define LEVEL 5

int main() {
  init_grid(1 << LEVEL);
  size (1); //dimension consistency

  mu1 = 1.;
  mu2 = 1.;

  run();
}

event init (t = 0) {
  solid (cs, fs, - y + 0.3);

  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);
  restriction ({cs, fs, cm, fm});

  fraction (f, 0.1 - x);
}

event final (i = 10) {
  output_ppm (cm, file = "cm.png");
}