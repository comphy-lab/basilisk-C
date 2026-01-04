#include "navier-stokes/centered.h"

int main() {
  origin (-0.5, -0.5);
  init_grid (8);
  DT = 1. [0,1];
  TOLERANCE = 1e-7;
  stokes = true;
  run();
}

event init (i = 0) {
  refine (level == 3 && sq(x) + sq(y) < sq(0.25));
  a[] = {0,-1.};
}

event check (i = 1) {
  output_cells (stdout);
#if 1
  foreach()
    fprintf (stderr, "%g %g %g %.6f %.6f\n", x, y, p[],
	     fabs(u.x[]), fabs(u.y[]));
#else
  foreach_face(y)
    fprintf (stderr, "%g %g %g\n", x, y, uf.y[]);
#endif
  foreach()
    assert (fabs (u.y[]) < 1e-6);
}
