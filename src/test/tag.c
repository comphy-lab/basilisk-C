#include "utils.h"
#include "tag.h"

int main()
{
  size (1.[0]);
  init_grid (64);
  X0 = Y0 = -0.5;
  scalar t[];
  foreach()
    t[] = sin(4*pi*x)*cos(4.*pi*y) > 0.5;
  int n = tag (t);
  foreach() {
    if (t[] == 0)
      t[] = nodata;
    else
      fprintf (qerr, "%g %g %g\n", x, y, t[]);
  }
  output_ppm (t, file = "t.ppm", min = 1, max = n, n = 256);
}
