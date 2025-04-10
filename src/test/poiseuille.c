#include "navier-stokes/centered.h"

int main() {
  size (1. [0]); // space and time are dimensionless
  DT = HUGE [0];
  
  origin (0, -0.5);
  stokes = true;
  TOLERANCE = 1e-5;
  for (N = 8; N <= 64; N *= 2)
    run(); 
}

u.t[top] = dirichlet(0);
u.t[bottom] = dirichlet(0);

u.n[left] = neumann(0);
p[left] = dirichlet(y + 0.5);
u.n[right] = neumann(0);
p[right] = dirichlet(y + 0.5);

scalar un[];

event init (t = 0) {
  // we also check for hydrostatic balance
  a[] = {1.,1.};
  mu[] = {1.,1.};
  foreach()
    un[] = u.x[];
}

event logfile (t += 0.1; i <= 100) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
}

event profile (t = end) {
  printf ("\n");
  foreach()
    printf ("%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
  scalar e[];
  foreach()
    e[] = u.x[] - 0.5*(0.25 - y*y);
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}
