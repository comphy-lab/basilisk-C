/**
# MUSCL TVD
I borrow the idea of bcg.h to create van-Leer TVD which will employed by other header file to facilitate MULES VOF scheme. 
*/

#undef SEPS
#define SEPS 1e-30

double vanleer_limiter(double r)
{
  return (r + fabs(r))/(1 + fabs(r));
}

foreach_dimension()
void vanleer_flux_x (scalar f,
		    face vector uf,
		    face vector flux)
{

  /**
  We first define limiter phi and delta f minus and plus. 
  */

  vector phi[];

  foreach()
  {
    /**
    Note the form:
    `df.x[] = f[1,0] - f[0,0];
    r = df.x[-1,0]/(df.x[0,0] + SEPS);`
    does not work well since it ignores the initial boundary which will set by default instead of calculating it following the rules here.
    */
    double r = (f[0,0] - f[-1,0])/(f[1,0] - f[0,0] + SEPS);
    phi.x[] = vanleer_limiter(r);
  }

  /**
  Accordingly we can obtain left and right value and eventually the flux. Note function sign can only return 1 (pos) and -1 (neg or 0)*/

  foreach_face(x) {
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

/**
 Update f here, since MULES requires flux on each face, the two step cannot merge. 
*/
foreach_dimension()
void TVD_update_x (scalar f, face vector flux, double dt)
{
  foreach()
  {
    f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
  }
}

/**
Advection below is used by tests for convenience.
*/

trace
void TVD_advection (scalar * tracers, face vector u, double dt, int i)
{
  scalar f;
  for (f in tracers) {
    face vector flux[];

    void (*faceflux[dimension]) (scalar, face vector, face vector);
    void (*update[dimension]) (scalar, face vector, double);

    int d = 0;

    foreach_dimension()
    {
      faceflux[d] = vanleer_flux_x;
      update[d++] = TVD_update_x;
    }

    for (d = 0; d < dimension; d++)
    {
      faceflux[(i + d) % dimension](f, u, flux);
      update[(i + d) % dimension](f, flux, dt);
    }
  }
}