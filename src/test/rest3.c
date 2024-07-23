// lake at rest with emerged island (same as rest1.c but for implicit
// Saint-Venant)
#include "saint-venant-implicit.h"

int main()
{
  origin (-0.5, -0.5);
  init_grid (32);
  DT = 1;
  TOLERANCE = 1e-12;
  dry = 1e-6;
  run();
}

event init (i = 0)
{
  refine (level == 5 && x < 0.1 && y < 0.1);

  double a = 0.3, b = 100.;
  foreach() {
    zb[] = a*exp(-b*(x*x + y*y));
    h[] = max (0.1 - zb[], 0.);
  }
}

event logfile (i = 1)
{
  norm n = normf (uf.x);
  fprintf (stderr, "# %.10f %.10f %.10f\n", n.avg, n.rms, n.max);
  printf ("x h zb uf.x a.x alpha.grad.p\n");
  foreach_face (x)
    printf ("%g %g %g %g %g %g\n",
	    x, h[], zb[], uf.x[], 
	    a.x[], alpha.x[]*(p[] - p[-1])/Delta);
}
