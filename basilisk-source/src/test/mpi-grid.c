#include "grid/multigrid.h"

int main()
{
  foreach_dimension()
    periodic (right);

  L0 = 2.*pi;
  init_grid (8);
  output_cells (stdout);

  vector u[];
  scalar s[];
  foreach()
    u.x[] = u.y[] = s[] = sin(x)*cos(y)*cos(z);

  foreach()
    foreach_neighbor(2) {
#if 0
      fprintf (stderr, "%g %g %g %g %g %g\n", x, y, z, s[],
	       s[] - sin(x)*cos(y)*cos(z), u.y[]);
#else
      assert (fabs (s[] - sin(x)*cos(y)*cos(z)) < 1e-14);
      assert (fabs (u.x[] - sin(x)*cos(y)*cos(z)) < 1e-14);
      assert (fabs (u.y[] - sin(x)*cos(y)*cos(z)) < 1e-14);
#endif
    }
}
