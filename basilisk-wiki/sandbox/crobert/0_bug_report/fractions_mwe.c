/**
# Fraction minimal working exemple

This file is a minimal working example to reproduce a small bug in fractions.h. 
If the Level-Set function used to define the fraction nullifies on the diagonal
of a cell, the resulting volume fraction of the cell is set at 0, whereas it should be 0.5.

A correction is suggested on the [forum](https://groups.google.com/forum/#!topic/basilisk-fr/N_TT5yxrSd0) and in the file [fractions_cr.h](fractions_cr.h).

If the flag SOLVED is set to 0, current [fractions.h](/src/fractions.h) is used and the 
bug can be seen in the plots at the end of the file. If SOLVED is set to 1, 
[fractions_cr.h](fractions_cr.h) is used and one can see that the problem is solved, 
at least in this particular case.
*/

#define SOLVED 0


#include "grid/multigrid.h"
#if SOLVED
  #include "fractions_cr.h"
#else
  #include "fractions.h"
#endif

/**
We set LEVEL to 1 to have a small number of cells and better see the bug. */

const int LEVEL = 2;

int main() {
  /** The orgin is set at the center of the domain. */
  origin (- 0.5, - 0.5);
  N = 1 << LEVEL;
  init_grid (N);

  /**
  Three fractions are defined:

  * $f_x$ corresponding to the Level Set function $x$;
  * $f_{x+y}$ corresponding to the Level Set function $x+y$;
  * $f_{-x-y}$ corresponding to the Level Set function $-x-y$.

  All three should correspond to an half space.*/

  scalar fx[], fxy[], f_xy[];

  fraction (fx, x);
  fraction (fxy, x+y);
  fraction (f_xy, -x-y);

  /**
  The cell grid and values of the fractions are saved. */

  FILE * fp = fopen ("cells", "w");
  output_cells (fp);
  fclose (fp);
      
  fp = fopen ("facets_fx", "w");
  output_facets(fx, fp);
  fclose (fp);
  
  fp = fopen ("facets_fxy", "w");
  output_facets(fxy, fp);
  fclose (fp);
  
  fp = fopen ("facets_f_xy", "w");
  output_facets(f_xy, fp);
  fclose (fp);
        
  fp = fopen ("fvalues", "w");
  foreach()
    fprintf(fp,"%g %g %g %g %g \n", x, y, fx[], fxy[], f_xy[]);
  fclose (fp);
}

/**

## Ouputs

Each value of the fractions resulting from the three Level Set functions
are written in its cell :

 ~~~gnuplot Reconstruction of $\phi(x, y) = x$
 set terminal @PNG enhanced size 640,640 font ",8"
 set output 'plot1.png'
 unset key 
 unset border
 unset tics
plot \
  'cells' w l lc rgb "#7F7F7F", \
  'facets_fx' u 1:2 w l, \
  'fvalues' u 1:2:3 with labels
 ~~~
 
 ~~~gnuplot Reconstruction of $\phi(x, y) = x+y$
 set output 'plot3.png'
plot \
  'cells' w l lc rgb "#7F7F7F", \
  'facets_fxy' u 1:2 w l, \
  'fvalues' u 1:2:4 with labels
 ~~~
  
 ~~~gnuplot Reconstruction of $\phi(x, y) = -x-y$
 set output 'plot4.png'
plot \
  'cells' w l lc rgb "#7F7F7F", \
  'facets_f_xy' u 1:2 w l, \
  'fvalues' u 1:2:5 with labels
 ~~~  

The diagonals initialisations are not correct
*/
