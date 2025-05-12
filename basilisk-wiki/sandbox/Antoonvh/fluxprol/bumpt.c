/**
# `saint-venant.h` on a quadtree

This is the mostly the `bump.c` code from the [tutorial](/Tutorial).

Changes include:  

* Periodic boundaries 
 */
#include "saint-venant.h"

#define LEVEL 8

int main() {
  periodic (left);     // This should not matter
  periodic (bottom);   
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  run();
}

event init (t = 0) {
  foreach()
    h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
}

event graphs (i++) {
  char fname[99];
  sprintf (fname, "%s", GRIDNAME);
  static FILE * fp = fopen (fname, "w");
  stats s = statsf (h);
  fprintf (fp, "%g %g %g\n", t, s.min, s.max);
}

event end (t = 4); 

event adapt (i++) {
  adapt_wavelet ({h}, (double[]){4e-3}, maxlevel = LEVEL);
}
