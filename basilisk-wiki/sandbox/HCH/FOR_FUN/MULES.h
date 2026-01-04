/**
 This is a file to facilitate MULES, given a VOF tracer and face velocity, the tracer will eventually updates with flux limiter $\lambda$. For each face we first need to obtain the lower flux. A simple up-wind is employed.
*/

#undef SEPS
#define SEPS 1e-30

void upwind(scalar f, face vector uf, face vector flux)
{
  foreach_face()
  {
    flux.x[] = uf.x[] > 0 ? (uf.x[]*f[-1, 0]):(uf.x[]*f[0, 0]);
  }
}

/**
Then we compute the flux correction directly.
*/

double MULES_vanleer_limiter(double r)
{
  return (r + fabs(r))/(1 + fabs(r));
}

void corr_flux (scalar f,
		    face vector uf,
		    face vector flux)
{

  vector phi[];

  foreach()
  {
    foreach_dimension()
    {
      double r = (f[0] - f[-1])/(f[1] - f[0] + SEPS);
      phi.x[] = MULES_vanleer_limiter(r);
    }
  }

  foreach_face() {
    double q;
    if (uf.x[] > 0.) {
      // Upstream left: reconstruct from f[i-1] to face i-1/2
      q =  0.5 * phi.x[-1] * (f[0] - f[-1]);
    } else {
      // Upstream right: reconstruct from f[i] to face i-1/2
      q =  -0.5 * phi.x[0] * (f[1] - f[0]);
    }
    flux.x[] = q * uf.x[];
  }
}

void limiter(scalar f,
	       	face vector flux_low,
	       	face vector flux_corr,
	       	face vector lambda,
	       	double dt,
		int cir = 10
		)
{
  /**
  Initialize lambda with value 1
  */
  foreach_face()
  {
    lambda.x[] = 1;
  }

  /**
  $P$ and $Q$ are static along the iteration, we shall compute them first
  */

  scalar P_p[], P_n[];
  scalar Q_p[], Q_n[];

  foreach()
  {
    P_p[] = P_n[] = Q_p[] = Q_n[] = 0;

    foreach_dimension()
    {
      if(flux_corr.x[] > 0.)
      {
        P_p[] += flux_corr.x[];
      }
      else
      {
        P_n[] -= flux_corr.x[];
      }

      if(flux_corr.x[1] > 0.)
      {
        P_n[] += flux_corr.x[1];
      }
      else
      {
        P_p[] -= flux_corr.x[1];
      }
    }

    double fmax = f[], fmin = f[];
    foreach_neighbor(1)
    {
      fmax = f[] > fmax ? f[] : fmax;
      fmin = f[] <= fmin ? f[] : fmin;
    }

    double low_flux = 0;
    foreach_dimension()
    {
      low_flux += flux_low.x[1] - flux_low.x[];
    }

    Q_p[] = (Delta*cm[])*(fmax - f[])/dt + low_flux;
    Q_n[] = (Delta*cm[])*(f[] - fmin)/dt - low_flux;
  }


  /**
  Then we can dive into the loop. First update the $(\sum \lambda C_f)^+$ and $(\sum \lambda C_f)^-$
  */
  scalar lambda_p[], lambda_n[];


  for (int i = 0; i < cir; i++)
  {
    scalar P_plambda[], P_nlambda[];

    foreach()
    {
      P_plambda[] = 0; P_nlambda[] = 0;

      foreach_dimension()
      {
        if(flux_corr.x[] > 0.)
        {
          P_plambda[] += lambda.x[]*flux_corr.x[];
        }
        else
        {
          P_nlambda[] -= lambda.x[]*flux_corr.x[];
        }

        if(flux_corr.x[1] > 0.)
        {
          P_nlambda[] += lambda.x[1]*flux_corr.x[1];
        }
        else
        {
          P_plambda[] -= lambda.x[1]*flux_corr.x[1];
        }
      }

      /**
      Now the cell-centered limiter $\lambda^+-$
      */
      double lambdap = 0, lambdan = 0;
      lambdap = (Q_p[] + P_nlambda[])/(P_p[]+SEPS);
      lambdan = (Q_n[] + P_plambda[])/(P_n[]+SEPS);

      lambda_p[] = max(min(lambdap,1), 0);
      lambda_n[] = max(min(lambdan,1), 0);
    }

    foreach_face()
    {
      lambda.x[] = flux_corr.x[] > 0? min(lambda_n[-1],lambda_p[]) 
         	    : min(lambda_n[],lambda_p[-1]);
    }
  }

}

trace
void MULES_advection (scalar * tracers, face vector u, double dt)
{
  scalar f;
  for (f in tracers) {

    face vector flux_corr[];
    face vector flux_low[];
    face vector lambda[];

    upwind (f, u, flux_low);
    corr_flux (f, u, flux_corr);
    limiter(f, flux_low, flux_corr, lambda, dt);

    foreach()
      foreach_dimension()
        f[] += dt*(flux_low.x[] - flux_low.x[1] +
			lambda.x[]*flux_corr.x[] - 
			lambda.x[1]*flux_corr.x[1])/(Delta*cm[]);
			//1*flux_corr.x[] - 
			//1*flux_corr.x[1])/(Delta*cm[]);
  }
}

void upwind_advection (scalar * tracers, face vector u, double dt)
{
  scalar f;
  for (f in tracers) {

    face vector flux[];

    upwind (f, u, flux);

    foreach()
      foreach_dimension()
        f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
  }
}
/**
The dimensional split is not employed here, consequently the result is not that ideal in multidimensional test, still wonder how OpenFOAM cope with it. Any special treatment?
*/