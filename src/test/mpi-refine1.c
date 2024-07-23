/* See ../figures/mpi-refine1.svg */

#include "refine_unbalanced.h"

int main (int argc, char * argv[])
{
#if 1
  init_grid (16);
#else
  init_grid(1);
  foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }
  tree->dirty = true;

  refine_unbalanced (level < 4, NULL);

  mpi_partitioning();
#endif
  refine_unbalanced (level < 5 && y < 0.315 && (1. - x) < 0.438 && y > 0.25 &&
		     (1. - x) > 0.375, NULL);
  unrefine (y < 0.25 && (1. - x) > 0.5);
  refine_unbalanced (level < 5 && y < 0.315 && (1. - x) < 0.5 && y > 0.25 &&
		     (1. - x) > 0.438, NULL);

  scalar s[];
  foreach()
    s[] = 0.;
  boundary ({s});
}
