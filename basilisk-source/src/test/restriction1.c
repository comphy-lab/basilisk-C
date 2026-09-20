/**
# Restriction on rectangular domains */

#include "grid/multigrid.h"

scalar s[];

int main()
{
#if dimension == 2
  dimensions (3, 1);
#elif dimension == 3
  dimensions (3, 2, 1);
#endif
  init_grid (16);
  foreach()
#if dimension == 2
    s[] = x*y;
#elif dimension == 3
    s[] = x + 10.*y + 100.*z;
#endif
  restriction ({s});
  foreach_level (2, serial)
#if dimension == 2
    fprintf (stderr, "%g %g %g\n", x, y, s[]);
#elif dimension == 3
    fprintf (stderr, "%g %g %g %g\n", x, y, z, s[]);
#endif
}
