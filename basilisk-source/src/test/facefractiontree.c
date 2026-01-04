/**
# Face fraction on trees

This test verifies that the face fraction calculation works well also on tree
grids. */

#include "fractions.h"

int main (void) {
  origin (-0.5, -0.5);
  init_grid (1 << 2);
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
  foreach_face()
    fprintf (stderr, "%g %g %.3g\n", x, y, sf.x[]);
  fp = fopen ("facets", "w");
  output_facets (c, fp);
  fclose (fp);
}

/**
~~~gnuplot Reconstructed face fractions
set size ratio -1
unset key 
unset border
unset tics
plot 'cells' w l, 'facets' w l, 'log' u 1:2:3 with labels
~~~
*/
