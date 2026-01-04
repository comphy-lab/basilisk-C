/**
# MUSCL no split TVD
This TVD scheme performs well in one-dimensional tests but fails in multidimensional cases due to the neglect of dimensional splitting.
*/

#undef SEPS
#define SEPS 1e-30

double vanleer_limiter_nosplit(double r)
{
  return (r + fabs(r))/(1 + fabs(r));
}

trace
void vanleer_flux (scalar f,
		    face vector uf,
		    face vector flux)
{

  vector phi[];

  foreach()
  {
    foreach_dimension()
    {
      double r = (f[0,0] - f[-1,0])/(f[1,0] - f[0,0] + SEPS);
      phi.x[] = vanleer_limiter_nosplit(r);
    }
  }

  foreach_face() {
    double q;
    if (uf.x[] > 0.) {
      // Upstream left: reconstruct from f[i-1] to face i-1/2
      q = f[-1,0] + 0.5 * phi.x[-1,0] * (f[0,0] - f[-1,0]);
    } else {
      // Upstream right: reconstruct from f[i] to face i-1/2
      q = f[0,0] - 0.5 * phi.x[0,0] * (f[1,0] - f[0,0]);
    }
    flux.x[] = q * uf.x[];
  }
}

trace
void TVD_nosplit_advection (scalar * tracers, face vector u, double dt)
{
  scalar f;
  for (f in tracers) {
    face vector flux[];
    vanleer_flux (f, u, flux);
    foreach()
      foreach_dimension()
        f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
  }
}