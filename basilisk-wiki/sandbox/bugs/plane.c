/**
# Slanted plane in bview

Planes may not be drawn correctly,

![A non unitary normal vector makes the cell slices be oriented out of the plane](plane/buggy.png)

![A unit normal makes some of the cell slices go missing](plane/buggy_2.png)

Setting a finite value for the plane *alpha* works around this problem.

![`alpha = 1e-4`](plane/no_bug.png)
 */

#include "grid/multigrid3D.h"
#include "view.h"

scalar s[];

int main() {
  X0 = Y0 = Z0 = -0.5;
  init_grid (N);
  coord normal = {2, 3, 1};
  box();
  squares ("s", n = normal);
  save ("buggy.png");
  // Normalizing does not fix this issue
  normalize (&normal);
  box();
  squares ("s", n = normal);
  save ("buggy_2.png");
  
  box();
  squares ("s", n = normal, alpha = 1e-4);
  save ("no_bug.png");
}
