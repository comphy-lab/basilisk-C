/**
# Adaption of cm & cs 

In case of adaption in a metric + embed problem the solid variables
must be at the head of the list of adapted variables.

The actual implementation in which the solid variables cs[] and fs[]
are defined in [embed.h](/src/embed.h) does not allow an easy solution
since the solid variables are not defined at the moment in which
adapt_wavelet() is coded. If the solid variables were defined as the
metric ones (at [common.h](/src/common.h)) it would be easy to code.

![cm[] adaption](axi-adapt/cm.png)
*/

#include "embed.h"
#include "axi.h"
#include "run.h"

#define LEVEL 5

scalar f[];

int main() {
  init_grid (1 << LEVEL);
  run();
}

event init(t=0) {
  solid (cs, fs, - y + 0.3);

  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);

  foreach()
    f[] = exp(-sq(x-0.5)/0.01);
}

#if 1
event adapt (i++) {
    adapt_wavelet ((scalar *){f}, (double[]){1e-3}, LEVEL + 2);
}
#else //Workaround
event adapt (i++) {
  adapt_wavelet ((scalar *){f}, (double[]){1e-3}, LEVEL + 2, list = {cs,fs,cm,fm});
}
#endif

event final (i=1) {
  output_ppm (cm, file = "cm.png", n = 300);
  return 1;
}

