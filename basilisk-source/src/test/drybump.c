#include "grid/multigrid1D.h"
#if EXPLICIT
# include "saint-venant.h"
#else
# include "saint-venant-implicit.h"
#endif

#define LEVEL 10

int main()
{
  origin (-0.5);
  init_grid (1 << LEVEL);
  DT = 1e-1;
  run();
}

event init (i = 0)
{
  double a = 1., l2 = 200.;
  foreach() {
    h[] = a*exp(- l2*(x*x));
    x -= -0.25;
    zb[] = a*exp(- l2*(x*x));
    h[] = max(h[] - zb[], 0.);
  }
}

event logfile (i++) {
  stats s = statsf (h);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum dt\n");
  fprintf (stderr, "%g %d %g %g %.8f %g\n", t, i, s.min, s.max, s.sum, dt);
  assert (s.min >= 0.);
}

event outputfile (t <= 0.6; t += 0.6/8) {
  static int nf = 0;
  printf ("file: eta-%d\n", nf++);
  foreach()
    printf ("%g %g %g\n", x, h[], zb[]);
  printf ("\n");
}
