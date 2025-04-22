/**
Checks that "half-mapped" fine/coarse cells do not cause trouble with
boundary conditions. */

#define TRASH 1
#include "grid/quadtree.h"

int main() {
  init_grid (N);
  
  scalar s[];  
  face vector v[];
  // mask (y >  L0/4. ? top : none);
  mask (sq(x - 0.5) + sq(y - 0.5) > sq(0.3) ? top : none);
  
  output_cells (stdout);
  
  for (int l = depth(); l > 0; l--)
    foreach_halo (restriction,l)
      fprintf (stderr, "%g %g %d\n", x, y, l);
  fflush (stderr);
  
  foreach()
    s[] = 1.;
  boundary ({s});

  foreach_face()
    v.x[] = 1.;
  boundary ((scalar *){v});

  restriction ((scalar *){v});
}
