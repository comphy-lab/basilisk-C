#include "mypoisson.h"

struct Viscosity {
  vector u;
  face vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
#if EMBED
  void (* embed_stress_flux) (Point, vector, vector, coord *, coord *);
#endif // EMBED
};

#if AXI
// fixme: RHO here not correct
# define lambda ((coord){1., 1. + dt/RHO*(mu.x[] + mu.x[1] + \
					  mu.y[] + mu.y[0,1])/2./sq(y)})
#else // not AXI
# if dimension == 1
#   define lambda ((coord){1.})
# elif dimension == 2
#   define lambda ((coord){1.,1.})
# elif dimension == 3
#   define lambda ((coord){1.,1.,1.})
#endif
#endif

// Temporary placement for tangential face gradients

#ifndef EMBED
#define face_avg_gradient_t1_x(a,i)				\
  ((a[1,i-1] + a[1,i] - a[-1,i-1] - a[-1,i])/(4.*Delta))
#define face_avg_gradient_t2_x(a,i) \
  ((a[1,0,i-1] + a[1,0,i] - a[-1,0,i-1] - a[-1,0,i])/(4.*Delta))

#define face_avg_gradient_t1_y(a,i) \
  ((a[i-1,1] + a[i,1] - a[i-1,-1] - a[i,-1])/(4.*Delta))
#define face_avg_gradient_t2_y(a,i) \
  ((a[0,1,i-1] + a[0,1,i] - a[0,-1,i-1] - a[0,-1,i])/(4.*Delta))

#define face_avg_gradient_t1_z(a,i) \
  ((a[i-1,0,1] + a[i,0,1] - a[i-1,0,-1] - a[i,0,-1])/(4.*Delta))
#define face_avg_gradient_t2_z(a,i) \
  ((a[0,i-1,1] + a[0,i,1] - a[i-1,0,-1] - a[0,i,-1])/(4.*Delta))
#endif // EMBED
  
// Note how the relaxation function uses "naive" gradients i.e. not
// the face_gradient_* macros.

static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

#if EMBED
  void (* embed_stress_flux) (Point, vector, vector,
			      coord *, coord *) = p->embed_stress_flux;
#endif // EMBED

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif
  
  foreach_level_or_leaf (l) {
    coord c = {0., 0., 0.}, d = {0., 0., 0.};
#if EMBED
    if (embed_stress_flux)
      embed_stress_flux (point, u, mu, &c, &d);
#endif // EMBED
    foreach_dimension() {
      w.x[] = (dt*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
               #if dimension > 1
		   + mu.y[0,1]*(u.x[0,1] +
				face_avg_gradient_t1_x (u.y, 1)*Delta)
		   - mu.y[]*(- u.x[0,-1] +
			     face_avg_gradient_t1_x (u.y, 0)*Delta)
               #endif
	       #if dimension > 2
		   + mu.z[0,0,1]*(u.x[0,0,1] +
				  face_avg_gradient_t2_x (u.z, 1)*Delta)
		   - mu.z[]*(- u.x[0,0,-1] +
			     face_avg_gradient_t2_x (u.z, 0)*Delta)
               #endif
		   ) + (r.x[] - dt*c.x)*sq(Delta))/
	(sq(Delta)*(rho[]*lambda.x + dt*d.x) + dt*(2.*mu.x[1] + 2.*mu.x[]
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

static double residual_viscosity (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
  double maxres = 0.;

#if EMBED
  void (* embed_stress_flux) (Point, vector, vector, coord *, coord *) = p->embed_stress_flux;
#endif
    
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  foreach_dimension() {
    face vector taux[];
    foreach_face(x)
      taux.x[] = 2.*mu.x[]*face_gradient_x (u.x, 0);
    #if dimension > 1
      foreach_face(y)
	taux.y[] = mu.y[]*(face_gradient_y (u.x, 0) + 
			   face_avg_gradient_t1_x (u.y, 0));
    #endif
    #if dimension > 2
      foreach_face(z)
	taux.z[] = mu.z[]*(face_gradient_z (u.x, 0) + 
			   face_avg_gradient_t2_x (u.z, 0));
    #endif
    boundary_flux ({taux});
    foreach (reduction(max:maxres)) {
      double a = 0.;
      coord c = {0., 0., 0.}, d = {0., 0., 0.};
#if EMBED
      if (embed_stress_flux)
	embed_stress_flux (point, u, mu, &c, &d);
#endif // EMBED
      foreach_dimension()
	a += taux.x[1] - taux.x[];
      res.x[] = r.x[] - rho[]*lambda.x*u.x[] + dt*(a/Delta - (c.x + d.x*u.x[]));
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
  boundary (resl);
#else
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    coord c = {0., 0., 0.}, d = {0., 0., 0.};
#if EMBED
  if (embed_stress_flux)
    embed_stress_flux (point, u, mu, &c, &d);
#endif // EMBED    
    foreach_dimension() {
      res.x[] = r.x[] - rho[]*lambda.x*u.x[] +
        dt*(2.*mu.x[1,0]*face_gradient_x (u.x, 1)
	    - 2.*mu.x[]*face_gradient_x (u.x, 0)
        #if dimension > 1
	    + mu.y[0,1]* (face_gradient_y (u.x, 1) +
			  face_avg_gradient_t1_x (u.y, 1))
	    - mu.y[]*(face_gradient_y (u.x, 0) +
		      face_avg_gradient_t1_x (u.y, 0))
	#endif
        #if dimension > 2
	    + mu.z[0,0,1]*(face_gradient_z (u.x, 1) +
			   face_avg_gradient_t2_x (u.z, 1))
	    - mu.z[]*(face_gradient_z (u.x, 0) +
		      face_avg_gradient_t2_x (u.z, 0))
	#endif
	    )/Delta - dt*(c.x + d.x*u.x[]);
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
#endif
  return maxres;
}

#undef lambda

double TOLERANCE_MU = 0.; // default to TOLERANCE

trace
mgstats viscosity (struct Viscosity p)
{
  vector u = p.u, r[];
  scalar rho = p.rho;  
  foreach()
    foreach_dimension()
      r.x[] = rho[]*u.x[];

  face vector mu = p.mu;
  restriction ({mu, rho});

#if EMBED
  p.embed_stress_flux = u.x.boundary[embed] != antisymmetry ? embed_stress_flux : NULL;
#endif // EMBED
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p, p.nrelax, p.res,
		   minlevel = 1, // fixme: because of root level
                                 // BGHOSTS = 2 bug on trees
		   tolerance = TOLERANCE_MU);
}
