/**
# Height functions

~~~gnuplot
set size ratio -1
unset key
set xrange [ -0.492536 : 0.0947108 ]
set yrange [ -0.459911 : 0.0620860 ]
plot 'out' w l, '< grep hx log' u ($1+$3*$4):2, \
     '< grep hy log' u 1:($2+$3*$4), \
     '< grep kappa log' u 1:2:3 w labels
~~~
*/

#include "fractions.h"
#include "curvature.h"
#include "utils.h"

scalar c[];

int main()
{
  init_grid (16);
  X0 = -0.5 - 0.125;
  Y0 = -0.5 - 0.125;
  Z0 = -0.5 - 0.125;
#if TREE
  c.refine = c.prolongation = fraction_refine;
  refine (level == 4 && fabs (x) < 0.375 && fabs (y) < 0.375);
  refine (level <= 5 && fabs (x) < 0.31 && fabs (y) < 0.31);
#endif
  fraction (c, - (0.2 - sqrt(sq(x+0.2) + sq(y+0.2))));
  output_cells (stdout);
#if dimension == 2
  output_facets (c, stdout);
#endif
  vector h[];
  heights (c, h);
  scalar kappa[];
  curvature (c, kappa);
#if dimension == 2
  foreach() {
    for (int i = -1; i <= 1; i++) {
      if (h.x[0,i] != nodata)
	fprintf (stderr, "%g %g %g %g hx\n",
		 x, y + i*Delta, height(h.x[0,i]), Delta);
      if (h.y[i] != nodata)
	fprintf (stderr, "%g %g %g %g hy\n",
		 x + i*Delta, y, height(h.y[i]), Delta);
    }
    if (kappa[] != nodata)
      fprintf (stderr, "%g %g %g kappa\n", x, y, kappa[]);
  }
#endif
  stats s = statsf (kappa);
  fprintf (stderr, "kappa min: %g avg: %g stddev: %g max: %g\n",
	   s.min, s.sum/s.volume, s.stddev, s.max);
#if 0
  FILE * fp = popen ("gfsview2D -s", "w");
  output_gfs (fp);
#endif
}
