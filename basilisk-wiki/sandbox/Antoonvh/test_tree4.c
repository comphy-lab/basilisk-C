/**
# 4 Ghosts 

~~~gnuplot
set size square
plot 'out' w l t 'Leaf cells', 'log' w labels t 'Prolongation cells'
~~~

What works:

* No crash
* The addional layers of ghosts cells are allocated
* These cells are marked as `prolongation` cells

What does not work

* No prolongation when parents are also prolongation cells
* The prolongation cell whos children are also prolngation cells are not marked as such. 
*/

#define dimension 2
#define GRIDNAME "Tree4"
#include "tree4.h"

void test_tree4_methods() {
  tree_methods();
}

scalar s[];
int main() {
  L0 = 5;
  X0 = Y0 = -L0/2;
  init_grid (1<<3);
  refine (sq(x) + sq(y) < sq(1) && level < 5);
  refine (x < 0 && x > -2 && level < 4);
  foreach() 
    s[] = x + y;
  boundary ({s});
  output_cells();
  foreach_cell()
    if (is_prolongation(cell))
      fprintf (stderr, "%g %g %.2g\n", x, y, s[]);
}
