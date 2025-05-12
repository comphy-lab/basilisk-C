/**
# Particles in a shallow-water solver

Small adjustments make it compatible:

![Particles in the bump example](bumpp/mov.mp4)

Note that the definitions are not compatible with the
`predictor-corrector` time loop, but is targeted at `run.h` loops.
 */
#include "saint-venant-implicit.h" 
#define u uf 
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

int LEVEL = 7;

Particles parts;

int main (int argc, char * argv[]) {
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  run();
}

event init (i = 0) {
  foreach()
    h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
  parts = init_tp_square (n = 5, l = 0.5);
}

event mov (t += 0.05) {
  scatter (parts);
  squares ("h", max = 0.5, min = 0.);
  save ("mov.mp4");
}

event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-3}, LEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}

event stop (t = 10);
