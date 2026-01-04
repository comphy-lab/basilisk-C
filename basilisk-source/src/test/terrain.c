#include "terrain.h"
#include "utils.h"

int main ()
{
  size (1.[0]);
  FILE * fp = popen ("xyz2kdt terrain", "w");
  for (double x = 0.; x <= 1.1; x += 0.005)
    for (double y = 0.; y <= 1.1; y += 0.005)
      fprintf (fp, "%g %g %g\n", x, y, sin(3.*pi*x)*cos(2.*pi*y));
  pclose (fp);

  for (int l = 4; l <= 7; l++) {
    init_grid (1 << l);
    scalar zb[];
    terrain (zb, "terrain", NULL);
    scalar e[];
    foreach()
      e[] = zb[] - sin(3.*pi*x)*cos(2.*pi*y);
    if (l == 7)
      output_field ({zb, e}, stdout, n = 128);
    norm n = normf (e);
    stats s = statsf (zb.nt);
    fprintf (stderr, "%d %.6f %.6f %.6f %g %g\n", l, n.avg, n.rms, n.max, 
	     s.min, s.max);
  }
}
