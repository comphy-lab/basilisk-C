#include "utils.h"

int main (int argc, char * argv[])
{
  init_grid (N);
  scalar s[];
  origin (-0.5, -0.5, -0.5);
  assert (restore (file = "restore-tree.dump", list = {s}));
  output_cells (stdout);
  double k = 1.;
  foreach()
    assert (s[] == sin(k*x)*cos(k*y));
}
