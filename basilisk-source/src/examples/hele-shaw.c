#include "hele-shaw.h"

#define MAXLEVEL 8

double mu1 = 2.5, mu2 = 1e-3, k = 0.1;

scalar f[];
scalar * tracers = {f};

p[left]  = dirichlet(1e-3);
p[right] = dirichlet(0);
f[left]  = 1.;

int main()
{
  size (7.5e-2);
  init_grid (1 << MAXLEVEL);
  run();
}

event init (i = 0) {
  foreach()
    f[] = (x < 1e-3)*(1. - 1e-3*(1. + noise()));
}

event logfile (i++)
{
  stats s = statsf (f);
  fprintf (stderr, "%d %g %d %g %g %g\n", 
	   i, t, mgp.i, s.sum, s.min, s.max);
}

event movies (t += 0.2)
{
  output_ppm (f, min = 0, max = 1, file = "f.mp4");

  output_ppm (u.y, file = "v.mp4");

  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, min = 0, max = MAXLEVEL, linear = false, file = "level.mp4");
}

event output (t += 10; t <= 50)
{
  static int nf = 0;
  printf ("file: field-%d\n", nf++);
  output_field ({p,f,u}, stdout, N);
}

event coefficients (i++)
{
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    beta.x[] = - k/(mu1 + clamp(ff,0,1)*(mu2 - mu1));
  }
}

#if TREE
event adapt (i++) {
  double tolerance = 0.02;
  adapt_wavelet ({f}, &tolerance, MAXLEVEL);
  event ("coefficients");
}
#endif
