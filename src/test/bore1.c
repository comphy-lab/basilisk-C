/**
# Simple Saint-Venant Riemann problem

This test is used by [Le Métayer et al,
2010](/src/references.bib#lemetayer2010) (section 6.1.1) to provide
the base solution for comparison with the Green-Naghdi equations.

## See also

* [Undular bores for the Green-Naghdi equations](bore.c) */

#include "grid/cartesian1D.h"
#include "saint-venant.h"

int main()
{
  X0 = -300.;
  L0 = 600.;
  G = 10.;
  N = 2048;
  run();
}

event init (i = 0)
{
  foreach()
    h[] = x < 0. ? 1.8 : 1.;
}

event output (t = 48) {
  foreach()
    fprintf (stdout, "%g %g %g\n", x, h[], u.x[]);
  fprintf (stdout, "\n");
}
