#include "grid/gpu/cartesian.h"
#include "saint-venant.h"

int LEVEL = 7;

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  origin (-0.5, -0.5);
  periodic (right);
  periodic (top);
  init_grid (1 << LEVEL);
  run();
}

#if _GPU && SHOW
event display (i++)
{
  output_ppm (h, fp = NULL, map = blue_white_red, n = 1024, spread = -1, fps = 30, linear = true);
}
#endif

event init (i = 0)
{
  double a = 1., b = 200.;
  foreach()
    h[] = 0.1 + a*exp(- b*(x*x + y*y));
}

#if !BENCHMARK
event logfile (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}
#endif

event outputfile (t <= 2.5; t += 2.5/8);
