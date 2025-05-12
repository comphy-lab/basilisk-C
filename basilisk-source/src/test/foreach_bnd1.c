/* Check that foreach_boundary() works with parallel multigrid */

#include "grid/multigrid.h"

int main() {
  init_grid (16);
  foreach_boundary (top)
    fprintf (qerr, "%g %g\n", x, y);
  foreach_boundary (bottom)
    fprintf (qerr, "%g %g\n", x, y);
  foreach_boundary (left)
    fprintf (qerr, "%g %g\n", x, y);
  foreach_boundary (right)
    fprintf (qerr, "%g %g\n", x, y);
}
