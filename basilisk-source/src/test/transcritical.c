#include "grid/cartesian1D.h"
#include "saint-venant.h"

int main()
{
  origin (-0.5, -0.5);
  init_grid (500);
  run();
}

u.n[left] = neumann(0);
u.n[right] = neumann(0);

event init (t = 0)
{
  double a = 0.25, k = pi/0.1;
  foreach() {
    zb[] = a*(cos(k*x) + 1.)*(fabs(x) < 0.1);
    u.x[] = 0.3;
    h[] = 1. - zb[];
  }
}

event logfile (t = 1.8) {
  foreach()
    fprintf (stderr, "%g %g %g %g\n", x, eta[], h[]*u.x[], zb[]);
}
