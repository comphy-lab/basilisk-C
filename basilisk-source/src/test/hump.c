// Small amplitude solitary wave interacting with a parabolic hump
#include "saint-venant.h"

#define MAXLEVEL 9
#define MINLEVEL 7

int main()
{
  init_grid (1 << MAXLEVEL);
  origin (-0.5, -1.);
  size (2.);
  run();
}

u.n[left]   = neumann(0);
u.n[right]  = neumann(0);
u.n[top]    = neumann(0);
u.n[bottom] = neumann(0);

event init (i = 0)
{
  foreach() {
    x += 0.5;
    zb[] = 0.8*exp(-5.*sq(x - 0.9) - 50.*y*y);
    h[] = (0.05 < x && x < 0.15 ? 1.01 : 1.) - zb[];
  }
}

event outputfile (t = {0.6, 0.9, 1.2, 1.5, 1.8})
{
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({eta}, stdout, N, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field ({l}, stdout, N);
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

event adapt (i++) {
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0.;

  astats s = adapt_wavelet ({eta}, (double[]){1e-4}, MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
