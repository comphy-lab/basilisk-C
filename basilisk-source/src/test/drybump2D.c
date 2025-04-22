#if IMPLICIT
# include "saint-venant-implicit.h"
#else
# include "saint-venant.h"
#endif

#define LEVEL 7

int main()
{
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  DT = 1e-1 [0,1];
  run();
}

#if IMPLICIT
q.n[right]  = neumann(0);
q.n[left]   = neumann(0);
q.n[top]    = neumann(0);
q.n[bottom] = neumann(0);
#else
u.n[right]  = neumann(0);
u.n[left]   = neumann(0);
u.n[top]    = neumann(0);
u.n[bottom] = neumann(0);
#endif

double a2 = 200.;

double terrain (double x, double y)
{
  x -= -0.25; y -= -0.25;
  double b = 0.5;
  return b*exp(-a2*(x*x + y*y));
}

void refine_zb (Point point, scalar zb)
{
  foreach_child()
    zb[] = terrain (x, y);
}

event init (i = 0)
{
  zb.refine = refine_zb; // updates terrain
  double b = 1.;
  foreach() {
    h[] = b*exp(-a2*(x*x + y*y));
    zb[] = terrain (x, y);
    h[] = max(h[] - zb[], 0.);
  }
}

event logfile (i++) {
  stats s = statsf (h);
#if IMPLICIT
  scalar u[];
  foreach()
    u[] = h[] > dry ? q.x[]/h[] : 0.;
  norm n = normf (u);
#else
  norm n = normf (u.x);
#endif
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %.8f %g %g %g\n", t, i, s.min, s.max, s.sum, 
	   n.rms, n.max, dt);
  //  assert (s.min >= 0.);
}

event outputfile (t <= 1.2; t += 1.2/8) {
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({h}, stdout, N, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field ({l}, stdout, N);
}

event adapt (i++) {
  // we do this so that wavelets use the default bilinear
  // interpolation this is less noisy than the linear + gradient
  // limiters used in Saint-Venant not sure whether this is better
  // though.
  scalar h1[];
  foreach()
    h1[] = h[];
  astats s = adapt_wavelet ({h1}, (double[]){1e-3}, LEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
