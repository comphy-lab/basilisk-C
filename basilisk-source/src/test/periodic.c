#include "poisson.h"
#include "utils.h"

int main()
{
  size (1.[0]);
  init_grid (N);
  scalar a[], b[], e[];

#if 0
  // this does not work with multigrid MPI or trees
  foreach_dimension() {
    a[right] = periodic();
  }
#else
  foreach_dimension()
    periodic (right);
#endif

  for (int n = 8; n <= 64; n *= 2) {
    init_grid (n);
    foreach() {
      a[] = 0.;
      b[] = 4.*dimension*sq(pi)*sin(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*z);
    }
    TOLERANCE = 1e-10;
    poisson (a, b);
    foreach()
      e[] = a[] + sin(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*z);
    stats s = statsf(e);
    foreach()
      e[] -= s.sum/s.volume;
    fprintf (stderr, "%d %g\n", n, statsf(e).max);
  }
  foreach()
    printf ("%g %g %g %g %g\n", x, y, z, a[], e[]);
}
