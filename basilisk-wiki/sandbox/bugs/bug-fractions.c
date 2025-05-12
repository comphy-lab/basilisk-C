/**
This bug could cause a problem when the fraction is calculated from a
distance function (level-set) by using "fractions()".  In this
example, an interface defined by the distance function crosses the
cell center of some cells.  The fraction of such cells is supposed to
be 0.5; however, "fractions()" gives 0 of the fraction to them. */

#include "fractions.h"

int main ()
{
  init_grid (16); // works with eight
  origin (0, 0, 0);
  L0 = 0.01;

  scalar levelset[], f[];
  
  foreach()
    levelset[] = x + y - 5e-3;
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = 0.25*(levelset[] + levelset[-1] + levelset[0,-1] + levelset[-1,-1]);
  fractions (phi, f);
  
  output_cells();
  foreach()
    fprintf (stderr, "%g %g %g\n", x, y, f[]);
}

/**
~~~gnuplot Volume fractions
set size ratio -1
unset key
unset xtics
unset ytics
set yrange [0:0.01]
plot 'out' w l, 'log' u 1:2:3 w labels, 5e-3 - x lw 3
~~~
*/
