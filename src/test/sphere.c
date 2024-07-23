/* similar to circle.c but in 3D (see also poisson3D.c) */

#include "grid/octree.h"
#include "utils.h"
#include "poisson.h"

double solution (double x, double y, double z)
{
  return sin(3.*pi*x)*sin(3.*pi*y)*sin(3.*pi*z);
}

int main()
{
  size (1.[0]);
  origin (-0.5, -0.5, -0.5);

  int depth = 6;
  init_grid(1);

  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y + z*z)));
  
  scalar a[], b[];

  foreach_dimension() {
    a[right] = dirichlet (solution(x, y, z));
    a[left]  = dirichlet (solution(x, y, z));
  }
  
  foreach() {
    a[] = 0.;
    b[] = -27.*pi*pi*sin(3.*pi*x)*sin(3.*pi*y)*sin(3.*pi*z);
  }

  poisson (a, b);
  
  double max = 0;
  foreach(reduction(max:max)) {
    double e = a[] - solution(x, y, z);
    if (fabs(e) > max) max = fabs(e);
    printf ("%g %g %g %g %g %g\n", x, y, z, a[], b[], e);
  }
  fprintf (stderr, "max error %g\n", max);
}
