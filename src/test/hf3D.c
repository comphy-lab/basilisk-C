#include "grid/octree.h"
#include "fractions.h"
#include "curvature.h"
#include "utils.h"

scalar f[];

int main()
{
  init_grid (16);
  origin (-0.5 - 0.125, -0.5 - 0.125, -0.5 - 0.125);
#if TREE
  refine (level == 4 && fabs (x) < 0.375 && fabs (y) < 0.375);
  refine (level <= 5 && fabs (x) < 0.31 && fabs (y) < 0.31);
  f.refine = f.prolongation = fraction_refine;
#endif
  fraction (f, - (0.2 - sqrt(sq(x+0.2) + sq(y+0.2) + sq(z))));
  vector h[];
  heights (f, h);
  scalar kappa[];
  curvature (f, kappa);
#if 0
  foreach() {
    if (h.x[] != nodata)
      fprintf (stderr, "%g %g %g %g hx\n", x, y, height(h.x[]), Delta);
    if (h.y[] != nodata)
      fprintf (stderr, "%g %g %g %g hy\n", x, y, height(h.y[]), Delta);
    if (kappa[] != nodata)
      fprintf (stderr, "%g %g %g kappa\n", x, y, kappa[]);
  }
#endif
  stats s = statsf (kappa);
  fprintf (stderr, "kappa min: %g avg: %g stddev: %g max: %g\n",
	   s.min, s.sum/s.volume, s.stddev, s.max);
#if 0
  FILE * fp = popen ("gfsview3D -s hf.gfv", "w");
  output_gfs (fp);
#endif
}
