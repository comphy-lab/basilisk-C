/**
# Fraction minimal working exemple 

This file is a minimal (not) working example to illustrate a small bug in
fractions.h. If the Level-Set function used to define the fraction is 0
exactly at a cell border, the two cells sharing this border may be set at 0,
whereas only one should.
[fractions_Q.h](/sandbox/qmagdelaine/geometry/fractions_Q.h) proposes a
correction. I am not sure my correction is good, see
[my message on the forum](https://groups.google.com/forum/#!topic/basilisk-fr/2mXqAAjiTdg)
for a more detailed explanation of the bug.

If the flag SOLVED is set to 0, current [fractions.h](/src/fractions.h) is
used and the bug can be seen in the plots at the end of the file. If SOLVED is
set to 1, [fractions_Q.h](/sandbox/qmagdelaine/geometry/fractions_Q.h) is used
and one can see that the problem is solved, at least in this particular case.
*/
#define SOLVED 0


#include "grid/multigrid.h"
#if SOLVED
  #include "fractions_Q.h"
#else
  #include "fractions.h"
#endif

/**
We set LEVEL to 2 to have a small number of cells and better see the bug. */

const int LEVEL = 2;

int main() {
  /** The orgin is set at the center of the domain. */
  origin (- 0.5, - 0.5);
  N = 1 << LEVEL;
  init_grid (N);

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

  FILE * fp = fopen ("cells", "w");
  output_cells (fp);
  fclose (fp);

  fp = fopen ("fvalues", "w");
  foreach()
    fprintf(fp,"%g %g %g %g %g %g \n", x, y, fx[], fy[], f_x[], f_y[]);
  fclose (fp);
}

/**

## Ouputs

Each value of the fractions resulting from the four Level Set functions
are written in its cell :

 ~~~gnuplot Reconstruction of $\phi(x, y) = x$
 set size 0.3
 unset key 
 unset border
 unset tics
 plot 'cells' w l, 'fvalues' u 1:2:3 with labels
 ~~~
 
 ~~~gnuplot Reconstruction of $\phi(x, y) = y$
 plot 'cells' w l, 'fvalues' u 1:2:4 with labels
 ~~~
 
 ~~~gnuplot Reconstruction of $\phi(x, y) = -x$
 plot 'cells' w l, 'fvalues' u 1:2:5 with labels
 ~~~
 
 ~~~gnuplot Reconstruction of $\phi(x, y) = -y$
 plot 'cells' w l, 'fvalues' u 1:2:6 with labels
 ~~~

All are correct except $\phi(x, y) = y$.
*/
