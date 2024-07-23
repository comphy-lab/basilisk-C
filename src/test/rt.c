/**
# Rayleigh--Taylor instability */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#if REDUCED
# include "reduced.h"
#endif

#define LEVEL 8

uf.n[left]   = 0;
p[left] = neumann(0.);
uf.n[right]  = 0;
p[right] = neumann(0.);

int main() {
  size (4);
  origin (-2, -2);
  init_grid (1 << LEVEL);
  mu1 = mu2 = 0.00313;
  rho1 = 1.225, rho2 = 0.1694;
#if REDUCED
  G.y = -9.81;
#else
  const face vector g[] = {0,-9.81};
  a = g;
#endif
  DT = 5e-3;
  TOLERANCE = 1e-6;
  run();
}

event init (t = 0) {
#if TREE  
  mask (x < -0.5 ? left : x > 0.5 ? right : none);
#endif
  fraction (f, 0.05*cos (2.*pi*x) + y);
}

event logfile (i++) {
  stats s = statsf (f);
  printf ("%g %d %g %g %g %g %d %d %d\n", 
	  t, i, dt, s.sum - 2., s.min, s.max - 1., mgp.i, mgpf.i, mgu.i);
  assert (s.min >= -1e-10 && s.max <= 1. + 1e-10);
  assert (fabs (s.sum - 2.) < 9e-7);
}

event interface (t = {0,0.2,0.4,0.8}) {
  output_facets (f, stderr);
  FILE * fp = fopen ("levels", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_field ({l}, fp);
  fclose (fp);
}

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen("gfsview2D -s rt.gfv", "w");
  output_gfs (fp);
}
#endif

#if TREE
event adapt (i++) {
  coord momb = {0};
  foreach()
    foreach_dimension()
      momb.x += dv()*rho(f[])*u.x[];
  
  adapt_wavelet ({f}, (double[]){5e-3}, LEVEL);

  coord moma = {0};
  foreach()
    foreach_dimension()
      moma.x += dv()*rho(f[])*u.x[];

  /**
  We check that each component of the momentum is conserved. */
  
  foreach_dimension()
    assert (fabs(momb.x - moma.x) < 1e-10);
}
#endif
