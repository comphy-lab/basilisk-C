/**
Check refinement of radial metric. */

#include "radial.h"
#include "run.h"

int main() {
  init_grid (4);
  run();
}

event init (i = 0) {
  refine (sq(x - 0.5) + sq(y - 0.5) < sq(0.3) && level < 3);
  foreach() {
    //    fprintf (stderr, "%g %g %g %g\n", x, y, cm[], r*dtheta/L0);
    assert (fabs(cm[] - r*dtheta/L0) < 1e-12);
  }
  foreach_face (x) {
    // fprintf (stderr, "%g %g %g %g\n", x, y, fm.x[], r*dtheta/L0);
    assert (fabs(fm.x[] - r*dtheta/L0) < 1e-12);
  }
}
