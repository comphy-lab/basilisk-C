/**
# Iterative convergence on a multigrid of is much worse than it is on a similar tree.

For this projection test case, the solver is unable to converge on a
multigrid, wheareas it works well on a quadtree.

Comment: the root of the problem here is that the right-hand-side is
non-zero on multigrids whereas it is indeed zero (as it should) on
quadtrees. So the problem seems to be that periodic boundary
conditions are not properly imposed for face vectors on multigrids... */

#include "grid/multigrid.h" // Comment for quadtree
#include "poisson.h"

int main() {
  init_grid (N);
  foreach_dimension()
    periodic (left);
  face vector v[];
  foreach_face()
    v.x[] = noise();
  scalar a[], div[];
  foreach() {
    a[] = div[] = 0;
    foreach_dimension()
      div[] = (v.x[1] - v.x[])/Delta;
  }
  poisson (a, div);
}
