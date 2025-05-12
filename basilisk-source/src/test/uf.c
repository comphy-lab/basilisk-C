/**
# Stability of the embedded face velocity interpolation

If a "naive" face interpolation is used (by the *face_value()* macro)
to compute face velocities in the [centered Navier--Stokes
solver](/src/navier-stokes/centered.h#viscous-term), an instability
can appear, due to the amplifications of velocity perturbations by the
third-order Dirichlet interpolation used to compute viscous fluxes.

This problem is avoided by using an embedded-fraction-weighted
interpolation of the face velocities (see [/src/embed.h]()). */

#include "embed.h"
#include "navier-stokes/centered.h"

int main()
{
  origin (-0.5, -0.5);
  periodic (right);
  periodic (top);
  
  stokes = true;
  DT = 2e-5;
  TOLERANCE = HUGE [*];
  NITERMIN = 10;
  N = 32;

  run();
}

event init (t = 0)
{
  double eps = L0/(1 << 7)/1000.;
  solid (cs, fs, union (y - L0/4. + eps, - L0/4 + eps - y));

  mu = fm;

  /**
  The boundary condition is zero velocity on the embedded boundary. */

  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);

  foreach()
    u.y[] = 1.;
}

event logfile (i++; i <= 100)
{
  fprintf (stderr, "%d %d %d %d %d %d %d %.3g %.3g %.3g %.3g\n",
	   i,
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resa*dt, mgu.resa, normf(u.y).max, normf(p).max);

  foreach()
    if (x < -L0/2. + L0/N)
      printf ("%g %g %g %g %g\n", y, u.y[], p[], g.y[], cs[]);
  printf ("\n");
}

event profile (t = end)
{
  p.nodump = false;
  dump();
  assert (normf(u.y).max < 1e-3);
}
