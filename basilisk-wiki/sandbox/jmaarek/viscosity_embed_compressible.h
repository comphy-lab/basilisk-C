#define COMPRESSIBLE 1

#include "poisson.h"

struct Viscosity {
  face vector mu;
  scalar rho;
  double dt;
  face vector lambdav;
  double (* embed_flux) (Point, scalar, vector, double *);
};

#if AXI
# define lambda ((coord){1., 1. + dt/((rho[]+SEPS)/(cs[] + SEPS))*(mu.x[]/(fs.x[] + SEPS) + mu.x[1]/(fs.x[1] + SEPS) + \
              mu.y[]/(fs.y[] + SEPS) + mu.y[0,1]/(fs.y[0,1] + SEPS))/2./sq(y+SEPS), 0})

#if COMPRESSIBLE
# define gamma ((coord) {1., 0, 0})

# define gamma1 ((coord) {0.,       \
      - (lambdav.x[]/(fs.x[] + SEPS) + lambdav.x[1]/(fs.x[1] + SEPS) + \
         lambdav.y[]/(fs.y[] + SEPS) + lambdav.y[0,1]/(fs.y[0,1] + SEPS))/4./sq(y+SEPS), 0})
#endif
#else // !AXI
# define lambda ((coord){1.,1.,1.})
#endif

/**
## Relaxation function

This function solves for the correction $d\boldsymbol{a}$ in step 3 of
the previously described algorithm. It is analogous to the relaxation
function written for the [Poisson-Helmholtz
equation](/src/poisson.h#application-to-the-poissonhelmholtz-equation),
with the only difference that the viscous diffusion equation is
vectorial and yields a system of coupled scalar equations, the number
of which is the dimension of the numerical simulation. This function
is passed as an argument to the [multigrid
cycle](/src/poisson.h#multigrid-cycle). */

static void relax_diffusion (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

  #if COMPRESSIBLE
  (const) face vector lambdav = p->lambdav;
#if AXI
  scalar ur = u.y;
#endif
#endif

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif

  double (* embed_flux) (Point, scalar, vector, double *) = p->embed_flux;
  foreach_level_or_leaf (l, nowarning) {

    foreach_dimension() {
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.x, mu, &c) : 0.;
      c *= cs[];
      d *= cs[];
      		
      w.x[] = ((r.x[] - dt*c)*sq(Delta) + dt*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]

#if COMPRESSIBLE
       + (lambdav.x[1]*u.x[1] + lambdav.x[]*u.x[-1]
               #if dimension > 1
          + lambdav.x[1,0]*((u.y[1,1] + u.y[0,1])/4 -
             (u.y[1,-1] + u.y[0,-1])/4.)
          - lambdav.x[]*((u.y[0,1] + u.y[-1,1])/4 -
             (u.y[0,-1] + u.y[-1,-1])/4.)
               #endif
               #if dimension > 2
          + lambdav.x[1,0,0]*((u.z[1,0,1] + u.z[0,0,1])/4 -
            (u.z[1,0,-1] + u.z[0,0,-1])/4.)
          - lambdav.x[]*((u.z[0,0,1] + u.z[-1,0,1])/4 -
             (u.z[0,0,-1] + u.z[-1,0,-1])/4.)
              #endif
#if AXI
          + (lambdav.x[1]*(ur[1] + gamma.x*ur[])  - lambdav.x[]*(ur[-1] + gamma.x*ur[]))/2.*Delta/y)*y
#else
                            )
#endif
#endif
#if dimension > 1
             + mu.y[0,1]*(u.x[0,1] +
		(u.y[1,0] + u.y[1,1])/4. -
		(u.y[-1,0] + u.y[-1,1])/4.)
             - mu.y[]*(- u.x[0,-1] +
		(u.y[1,-1] + u.y[1,0])/4. -
		(u.y[-1,-1] + u.y[-1,0])/4.)
#endif
#if dimension > 2
             + mu.z[0,0,1]*(u.x[0,0,1] +
		(u.z[1,0,0] + u.z[1,0,1])/4. -
		u.z[-1,0,0] + u.z[-1,0,1])/4.)
             - mu.z[]*(- u.x[0,0,-1] +
		(u.z[1,0,-1] + u.z[1,0,0])/4. -
		(u.z[-1,0,-1] + u.z[-1,0,0])/4.)
#endif
             ))/
(sq(Delta)*(rho[]*lambda.x + dt*d) + dt*(2.*mu.x[1] + 2.*mu.x[]
                    #if dimension > 1
				   + mu.y[0,1] + mu.y[]
                    #endif
	            #if dimension > 2
				   + mu.z[0,0,1] + mu.z[]
	            #endif
				   ) + SEPS);
    }
  }
  
#if JACOBI
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = (u.x[] + 2.*w.x[])/3.;
#endif

#if TRASH
  vector u1[];
  foreach_level_or_leaf (l)
    foreach_dimension()
      u1.x[] = u.x[];
  trash ({u});
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = u1.x[];
#endif
}

static double residual_diffusion (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  double (* embed_flux) (Point, scalar, vector, double *) = p->embed_flux;
  vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
  double maxres = 0.;

#if COMPRESSIBLE
  (const) face vector lambdav = p->lambdav;
#if AXI
  scalar ur = u.y;
#endif
#endif

#if TREE
  /* conservative coarse/fine discretisation (2nd order) */

  /**
  We manually apply boundary conditions, so that all components are
  treated simultaneously. Otherwise (automatic) BCs would be applied
  component by component before each foreach_face() loop. */

  boundary ({u});

  foreach_dimension() {
    face vector taux[];
#if COMPRESSIBLE
    face vector tauc[];
#if AXI
    face vector axic[];
#endif
#endif
    foreach_face(x){
#if COMPRESSIBLE
#if AXI
      axic.x[] = lambdav.x[]*(ur[]+ur[-1])/2.;

#endif
      tauc.x[] = lambdav.x[]*(u.x[] - u.x[-1]
    #if dimension > 1
            + (u.y[0,1] + u.y[-1,1])/4 
            - (u.y[0,-1] + u.y[-1,-1])/4.
    #endif
    #if dimension > 2
            + (u.z[0,0,1] + u.z[-1,0,1])/4 
            - (u.z[0,0,-1] + u.z[-1,0,-1])/4.
    #endif
            )/Delta;
#endif
      taux.x[] = 2.*mu.x[]*(u.x[] - u.x[-1])/Delta;
    }
    #if dimension > 1
      foreach_face(y){
	taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] + 
			   (u.y[1,-1] + u.y[1,0])/4. -
			   (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;}
    #endif
    #if dimension > 2
      foreach_face(z){
	taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] + 
			   (u.z[1,0,-1] + u.z[1,0,0])/4. -
			   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/Delta;}
    #endif

    //fprintf(ferr, "test\n");
    foreach (reduction(max:maxres), nowarning) {
      double a = 0.;
      foreach_dimension()
		  a += taux.x[1] - taux.x[];
      //fprintf(ferr, "%g %g %g\n", r.x[], rho[]*lambda.x*u.x[], dt*a/Delta);
      res.x[] = r.x[] - rho[]*lambda.x*u.x[] + dt*a/Delta;
      #if COMPRESSIBLE
#if AXI
      res.x[] += dt*((tauc.x[1] - tauc.x[] + 
          (axic.x[1] - axic.x[])/(y))/Delta
         + gamma1.x*ur[])*(y);
#else
      res.x[] += dt*(tauc.x[1]-tauc.x[])/Delta;
#endif
#endif
      if (embed_flux) {
		  double c, d = embed_flux (point, u.x, mu, &c);
	      res.x[] -= dt*(c + d*u.x[]);
      }
      if (fabs (res.x[]) > maxres)
	      maxres = fabs (res.x[]);
    }
    //fprintf(ferr, "test2\n");
  }
  boundary (resl);
#endif
  return maxres;
}

#if COMPRESSIBLE & AXI
#undef gamma
#undef gamma1
#endif
#undef lambda

double TOLERANCE_MU = 0.; // default to TOLERANCE

trace
mgstats viscosity (vector u, face vector mu, scalar rho, double dt,
		   int nrelax = 4, scalar * res = NULL, face vector lambdav)
{
  vector r[];
  foreach()
    foreach_dimension()
      r.x[] = rho[]*u.x[];

  /**
  We need $\mu$ and $\rho$ on all levels of the grid. */

#if COMPRESSIBLE
  restriction ({mu, lambdav, rho});
  struct Viscosity p = { mu, rho, dt, lambdav};
#else
  restriction ({mu,rho});
  struct Viscosity p = { mu, rho, dt};
#endif

  p.embed_flux = u.x.boundary[embed] != antisymmetry ? embed_flux : NULL;
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_diffusion, relax_diffusion, &p, nrelax, res,
		   minlevel = 1, // fixme: because of root level
                                  // BGHOSTS = 2 bug on trees
		   tolerance = TOLERANCE_MU ? TOLERANCE_MU : TOLERANCE);
}
