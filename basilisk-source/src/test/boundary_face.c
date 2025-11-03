/**
# Periodic boundary conditions on face fields

Checks on different grids in Makefile. */

#include "utils.h"

int main() {
  init_grid (2);
  foreach_dimension()
    periodic (left);
  face vector v[];
  foreach_face()
    v.x[] = (x + 1)*(3.*y + 7.);
  scalar a[], div[];
  foreach() {
    a[] = div[] = 0;
    foreach_dimension()
      div[] = (v.x[1] - v.x[])/Delta;
    fprintf (stderr, "c %g %g %g\n", x, y, div[]);
    fprintf (stderr, "F %g %g %g\n", x + Delta/2., y, v.x[1]);
    fprintf (stderr, "F %g %g %g\n", x - Delta/2., y, v.x[0]);
    fprintf (stderr, "F %g %g %g\n", x, y + Delta/2., v.y[0,1]);
    fprintf (stderr, "F %g %g %g\n", x, y - Delta/2., v.y[0]);
  }
  output_cells();
}
