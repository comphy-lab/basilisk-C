/**
~~~gnuplot Cells and field
unset key
set size ratio -1
plot 'out' u 1:2 w l, 'log' u 1:2:3 w labels
~~~
*/

#include "grid/quadtree.h"

int main()
{
  L0 = 16;
  X0 = Y0 = 0;
  // coordinates of lower-left corner
  origin (-0.5, -0.5);
  init_grid (4); // Initialize a 4 x 4 grid
  scalar h[];
  foreach()
    h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
  boundary ({h});
  refine ((x > 12) && (y > 12) && (level < 3)); // Refine to top right corner
  unrefine ((x < 8) && (y < 8) && level >= 1); // Coarsen the bottom left corner

  output_cells();
  foreach()
    fprintf (stderr, "%g %g %g\n", x, y, h[]);
  
  free_grid();
}
