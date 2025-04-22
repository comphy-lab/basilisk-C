/**
# Fractions in marginal cases

This checks that straight interfaces falling exactly on grid lines
give the correct fractions. */

#include "fractions.h"

int main() {
  
  /**
  The origin is set at the center of the domain. */
  
  origin (- 0.5, - 0.5);
  init_grid (4);

  /**
  Four fractions are defined:

  * $f_x$ corresponding to the Level Set function $x$;
  * $f_y$ corresponding to the Level Set function $y$;
  * $f_{-x}$ corresponding to the Level Set function $-x$;
  * $f_{-y}$ corresponding to the Level Set function $-y$.

  All four should correspond to an half space.*/

  scalar fx[], f_x[], fy[], f_y[];

  fraction (fx, x);
  fraction (fy, y);
  fraction (f_x, -x);
  fraction (f_y, -y);
  
  /**
  The cell grid and values of the fractions are saved. */

  output_cells (stdout);

  foreach()
    fprintf (stderr, "%g %g %g %g %g %g \n", x, y, fx[], fy[], f_x[], f_y[]);
}

/**
## Ouputs

 ~~~gnuplot Reconstruction of $\phi(x, y) = x$
 set term @SVG size 100,100
 unset key 
 unset border
 unset tics
 plot 'out' w l, 'log' u 1:2:3 with labels
 ~~~
 
 ~~~gnuplot Reconstruction of $\phi(x, y) = y$
 set term @SVG size 100,100
 unset key
 plot 'out' w l, 'log' u 1:2:4 with labels
 ~~~
 
 ~~~gnuplot Reconstruction of $\phi(x, y) = -x$
 set term @SVG size 100,100
 unset key
 plot 'out' w l, 'log' u 1:2:5 with labels
 ~~~
 
 ~~~gnuplot Reconstruction of $\phi(x, y) = -y$
 set term @SVG size 100,100
 unset key
 plot 'out' w l, 'log' u 1:2:6 with labels
 ~~~
*/
