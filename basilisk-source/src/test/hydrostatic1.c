// similar to hydrostatic.c but with a linear density profile
// (i.e. a quadratic pressure profile).

#include "navier-stokes/centered.h"

int main() {
  origin (-0.5, -0.5);
  init_grid (8);
  DT = 1.;
  TOLERANCE = 1e-6;
  stokes = true;
  run();
}

p[top] = dirichlet(0);

event init (i = 0) {
  refine (level == 3 && sq(x) + sq(y) < sq(0.25));
  a[] = {0,-1.};
  alpha = new face vector;
  foreach_face()
    alpha.x[] = 1./(0.51 - y);
}

event check (i = 1) {
  output_cells (stdout);
#if 1
  foreach()
    fprintf (stderr, "%g %g %g %.6f %.6f\n", x, y, p[], u.x[], u.y[]);
#else
  foreach_face(y)
    fprintf (stderr, "%g %g %g\n", x, y, uf.y[]);
#endif
}
