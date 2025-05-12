#include "poisson.h"

struct Viscosity {
  vector u;
  face vector mu;
  scalar rho;
  double dt;
  int nrelax;
#if COMPRESSIBLE
  face vector lambdav;
#endif
  scalar * res;
};

#if AXI
# define lambda ((coord){1., 1. + dt/rho[]*(mu.x[] + mu.x[1] + \
					    mu.y[] + mu.y[0,1])/2./sq(y)})
#if COMPRESSIBLE

# define gamma ((coord) {1., 0})

# define gamma1 ((coord) {0.,				\
			- (lambdav.x[] + lambdav.x[1] +	\
			   lambdav.y[] + lambdav.y[0,1])/4./sq(y)})
#endif
#else // not AXI
# if dimension == 1
#   define lambda ((coord){1.})
# elif dimension == 2
#   define lambda ((coord){1.,1.})
# elif dimension == 3
#   define lambda ((coord){1.,1.,1.})
#endif
#endif

static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
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
  foreach_level_or_leaf (l) {
    foreach_dimension() {
      w.x[] = (dt/rho[]*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
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
			    + (lambdav.x[1]*(ur[1] + gamma.x*ur[])  - 
			       lambdav.x[]*(ur[-1] + gamma.x*ur[]))/2.*Delta/y)*y
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
					  (u.z[-1,0,0] + u.z[-1,0,1])/4.)
			   - mu.z[]*(- u.x[0,0,-1] +
				     (u.z[1,0,-1] + u.z[1,0,0])/4. -
				     (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
               #endif
			   ) + r.x[]*sq(Delta))/
	(sq(Delta)*lambda.x + dt/rho[]*(2.*mu.x[1] + 2.*mu.x[]
                                   #if COMPRESSIBLE
					+ (lambdav.x[1] + lambdav.x[] 
                                       #if AXI
					   - (lambdav.x[1] - lambdav.x[])*Delta/2/y*gamma.y 
					   -  gamma1.x*sq(Delta))*y
                                       #else
			                  )
				       #endif   
                                    #endif
                                    #if dimension > 1
					+ mu.y[0,1] + mu.y[]
                                    #endif
			            #if dimension > 2
				      + mu.z[0,0,1] + mu.z[]
			            #endif
			     ));
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

static double residual_viscosity (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
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
  foreach_dimension() {
    face vector taux[];
#if COMPRESSIBLE
    face vector tauc[];
#if AXI
    face vector axic[];
#endif
#endif
    foreach_face(x) {
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

      taux.x[] =  2.*mu.x[]*(u.x[] - u.x[-1])/Delta;
    }
    #if dimension > 1
      foreach_face(y)
	taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] + 
			   (u.y[1,-1] + u.y[1,0])/4. -
			   (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;
    #endif
    #if dimension > 2
      foreach_face(z)
	taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] + 
			   (u.z[1,0,-1] + u.z[1,0,0])/4. -
			   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/Delta;
    #endif
    boundary_flux ({taux});
    foreach (reduction(max:maxres)) {
      double d = 0.;
      foreach_dimension()
	d += taux.x[1] - taux.x[];
      res.x[] = r.x[] - lambda.x*u.x[] + dt/rho[]*d/Delta;
#if COMPRESSIBLE
#if AXI
      res.x[] += dt/rho[]*((tauc.x[1] - tauc.x[] + 
			    (axic.x[1] - axic.x[])/y)/Delta
			   + gamma1.x*ur[])*y;
#else
      res.x[] += dt/rho[]*(tauc.x[1]-tauc.x[])/Delta;
#endif
#endif
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
  boundary (resl);
#else
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres))
    foreach_dimension() {
      res.x[] = r.x[] - lambda.x*u.x[] +
        dt/rho[]*(2.*mu.x[1,0]*(u.x[1] - u.x[])
		  - 2.*mu.x[]*(u.x[] - u.x[-1])
#if COMPRESSIBLE
		  + (lambdav.x[1]*(u.x[1] - u.x[])
		  - lambdav.x[]*(u.x[] - u.x[-1])
    #if dimension > 1
		  + lambdav.x[1]*((u.y[1,1] + u.y[0,1])/4 -
				  (u.y[1,-1] + u.y[0,-1])/4.)
		  - lambdav.x[]*((u.y[0,1] + u.y[-1,1])/4 -
				 (u.y[0,-1] + u.y[-1,-1])/4.)
    #endif
    #if dimension > 2
		  + lambdav.x[1]*((u.z[1,0,1] + u.z[0,0,1])/4 -
				  (u.z[1,0,-1] + u.z[0,0,-1])/4.)
		  - lambdav.x[]*((u.z[0,0,1] + u.z[-1,0,1])/4 -
				 (u.z[0,0,-1] + u.z[-1,0,-1])/4.)
   #endif
        #if AXI
		     + (lambdav.x[1]*(ur[1] + ur[])/2.
			- lambdav.x[]*(ur[-1] + ur[])/2.)*Delta/y
		     + gamma1.x*ur[]*sq(Delta))*y
        #else
		     )
        #endif
#endif
        #if dimension > 1
		  + mu.y[0,1]*(u.x[0,1] - u.x[] +
			       (u.y[1,0] + u.y[1,1])/4. -
			       (u.y[-1,0] + u.y[-1,1])/4.)
		  - mu.y[]*(u.x[] - u.x[0,-1] +
			    (u.y[1,-1] + u.y[1,0])/4. -
			    (u.y[-1,-1] + u.y[-1,0])/4.)
	#endif
        #if dimension > 2
		  + mu.z[0,0,1]*(u.x[0,0,1] - u.x[] +
				 (u.z[1,0,0] + u.z[1,0,1])/4. -
				 (u.z[-1,0,0] + u.z[-1,0,1])/4.)
		  - mu.z[]*(u.x[] - u.x[0,0,-1] +
			    (u.z[1,0,-1] + u.z[1,0,0])/4. -
			    (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
	#endif
		  )/sq(Delta);
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
#endif
 return maxres;
}

#if COMPRESSIBLE & AXI
#undef gamma
#undef gamma1
#endif
#undef lambda

mgstats viscosity (struct Viscosity p)
{
  vector u = p.u, r[];
  foreach()
    foreach_dimension()
      r.x[] = u.x[];

  face vector mu = p.mu;
  scalar rho = p.rho;
#if COMPRESSIBLE
  face vector lambdav = p.lambdav;
  restriction ({mu, lambdav, rho});
#else
  restriction ({mu,rho});
#endif
  
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p, p.nrelax, p.res);
}

mgstats viscosity_explicit (struct Viscosity p)
{
  vector u = p.u, r[];
  mgstats mg = {0};
  mg.resb = residual_viscosity ((scalar *){u}, (scalar *){u}, (scalar *){r}, &p);
  foreach()
    foreach_dimension()
      u.x[] += r.x[];
  boundary ((scalar *){u});
  return mg;
}
