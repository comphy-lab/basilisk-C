/**
# Restriction on rectangular domains */

#include "grid/multigrid.h"

scalar s[];

int main()
{
  dimensions (3, 1);
  init_grid (16);
  foreach()
    s[] = x*y;
  restriction ({s});
  foreach_level (2, serial)
    fprintf (stderr, "%g %g %g\n", x, y, s[]);
}
