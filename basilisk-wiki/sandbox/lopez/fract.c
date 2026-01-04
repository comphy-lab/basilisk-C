/**
# Test for reconstructed face fractions. */

//#include "grid/cartesian.h"
#include "fractions.h"
#include "src/fracface.h"

int main()
{
  origin (-0.5, -0.5);
  init_grid (4);
#if TREE
  refine (level == 2 && x < 0.25);
#endif
  FILE * fp = fopen ("cells", "w");
  output_cells (fp);
  fclose (fp);
 
  vertex scalar phi[];
  foreach_vertex()
    phi[] = 0.1 - sq(x) - sq(y);
  
  scalar c[];
#if TREE
  c.prolongation = fraction_refine;
#endif
  face vector sf[];
  fractions (phi, c);
  face_fraction (c, sf);
  fp = fopen ("log", "w");
  foreach_face()
    fprintf(fp,"%g %g %g \n", x, y, sf.x[]);
  fclose (fp);
  fp = fopen ("facet", "w");
  output_facets (c, fp);
  fclose (fp);
}

/**
~~~gnuplot Reconstructed face fractions
set size ratio -1
unset key 
unset border
unset tics
plot 'cells' w l, 'facet' w l, 'log' u 1:2:3 with labels
~~~
*/
