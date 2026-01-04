#include "navier-stokes/centered.h"

double tEnd = 1.;

int main() {
  run();
}

event logfile (i++) {
  fprintf (stderr, "%d %g %g\n", i, t, dt);
  assert (dt > 0.01);
}

event every30 (t += tEnd/30.; t <= tEnd) {
  fprintf (stderr, "every30 %g\n", t);
}

event end (t = tEnd) {
  fprintf (stderr, "tEnd %g\n", t);  
}
