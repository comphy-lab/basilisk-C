#include "advection.h"
#include "vof.h"

scalar c[];
scalar * interfaces = {c}, * tracers = NULL;
int MAXLEVEL;

int main()
{
  // coordinates of lower-left corner
  origin (-0.5, -0.5);
  for (MAXLEVEL = 5; MAXLEVEL <= 7; MAXLEVEL++) {
    init_grid (1 << MAXLEVEL);
    run ();
  }
}

#define circle(x,y) (sq(0.1) - (sq(x-0.25) + sq(y)))

event init (i = 0)
{
  fraction (c, circle(x,y));
}

#define end 0.785398

event velocity (i++) {
#if TREE
  double cmax = 1e-2;
  adapt_wavelet ({c}, &cmax, MAXLEVEL, list = {c});
#endif

  double a = -8.;
  trash ({u});
  foreach_face(x) u.x[] = - a*y;
  foreach_face(y) u.y[] =   a*x;
}

event logfile (t = {0,end}) {
  stats s = statsf (c);
  fprintf (stderr, "# %f %.12f %f %g\n", t, s.sum, s.min, s.max);
}

event interface (t += end/10.) {
  static FILE * fp = fopen ("interface", "w");
  if (N == 64) {
    output_facets (c, fp);
    if (t == end)
      output_cells (fp);
  }
}

event field (t = end) {
  scalar e[];
  fraction (e, circle(x,y));
  foreach()
    e[] -= c[];
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
  if (N == 64)
    output_field ({e}, stdout, N, linear = false);
}
