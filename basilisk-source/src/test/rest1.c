// lake at rest with emerged island
#include "saint-venant.h"

int main()
{
  origin (-0.5, -0.5);
  init_grid (32);
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
  norm n = normf (u.x);
  fprintf (stderr, "# %.10f %.10f %.10f\n", n.avg, n.rms, n.max);
  printf ("x y h zb u.x u.y eta\n");
  foreach ()
    printf ("%g %g %g %g %.3g %.3g %.3g\n", x, y, h[], zb[], 
	    u.x[] < 1e-10 ? 0. : u.x[], 
	    u.y[] < 1e-10 ? 0. : u.y[], 
	    h[] > dry ? h[] + zb[] - 0.1 : nodata);
}
