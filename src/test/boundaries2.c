/**
# Boundary conditions with 2 ghost layers */

#define BGHOSTS 2
#include "grid/multigrid.h"

int main()
{
  init_grid (4);
  output_cells (stdout);

  vector v[];
  scalar s[];
  face vector f[];

  s[left] = dirichlet(x + y);
  s[top] = dirichlet(x + y);
  s[right] = dirichlet(x + y);
  s[bottom] = dirichlet(x + y);

  v.t[bottom] = dirichlet(0);
  f.n[bottom] = 0.;

  foreach() {
    s[] = x + y;
    v.x[] = v.y[] = 1.;
  }
  foreach_face()
    f.x[] = 1.;

  foreach (serial)
    foreach_neighbor()
      fprintf (stderr, "%g %g %g %g %g\n", x, y, s[], v.x[], v.y[]);
  foreach_face (x, serial)
    for (int i = -2; i <= 2; i++)
      fprintf (stderr, "%g %g %g\n", x, y + i*Delta, f.x[0,i]);
  foreach_face (y, serial)
    for (int i = -2; i <= 2; i++)
      fprintf (stderr, "%g %g %g\n", x + i*Delta, y, f.y[i]);
}
