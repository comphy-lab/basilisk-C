/**
# Distance to a volume fraction field

`distance_to_surface` can be used to obtain a distance field.

![For example, to compute the distance to the green surface.](example_frac-dist/image.png)
*/

#include "grid/octree.h"
#include "frac-dist.h"
#include "view.h"

scalar cs[];
face vector fs[];

int main() {
  cs.prolongation = fraction_refine;
  periodic (left);
  periodic (right);
  L0 = 2*pi;
  X0 = Z0 = -L0/2;
  init_grid (8);
  
  vertex scalar phi[], dv[];
  scalar d[];
  do {
    foreach_vertex()
      phi[] = y + 0.4*sin(x) * cos(2*z) - 1.1;
    fractions (phi, cs, fs);
    boundary ({cs});
    distance_to_surface (cs, fs, d, dv);
  } while (adapt_wavelet ({d}, (double[]){0.005}, 8, 5).nf);
  
  view (ty = -0.3, theta = 0.4, phi = 0.4);
  draw_vof ("cs", "fs", fc = {0.2, 0.8, 0.2});
  isosurface ("d", 1);
  squares ("dv", alpha = -L0/2);
  cells(alpha = -0.5);
  save ("image.png");
}
