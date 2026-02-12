/**
This file is a modification of viscosity-embed.h without the supposition that $\mu$ should be constant. It implements this new term with an explicit version */


/**
# Viscous solver with embedded boundaries

We consider the Stokes equations
$$
\rho \mathbf{u}_t = \rho \mathbf{g} +
\nabla \cdot [\mu ( \nabla \mathbf{u} + \nabla^T \mathbf{u})]
$$
with $\mathbf{g}$ the acceleration and $T$ the transpose.

We solve implictly this equation using the multigrid solver
$$
dt\nabla \cdot [\mu ( \nabla \mathbf{u}^{n+1} + \nabla^T \mathbf{u}^{n+1})] - (\rho + \mathbf{\lambda}) \mathbf{u}^{n+1}
+ \underbrace{(\rho  \mathbf{u}^{n+1}+ \mathbf{g}\, dt) }_{r_i} = \mathbf{0}
$$
$\mathbf{\lambda}$ is a possible extra term due to the metric. 

The flux coming from the solid is obtained like the viscous force in embed.h*/

#include "poisson_perso.h"

struct Viscosity {
  face vector mu;
  scalar rho;
  double dt;
#if EMBED
  double (* embed_flux) (Point, scalar, vector, double *);
#endif // EMBED
};

#if AXI
// fixme: RHO here not correct
# define lambda ((coord){0.,dt*(mu.x[] + mu.x[1] +	\
					  mu.y[] + mu.y[0,1])/2./sq(y),0.})
#else // not AXI
# if dimension == 1
#   define lambda ((coord){0.})
# elif dimension == 2
#   define lambda ((coord){0.,0.})
# elif dimension == 3
#   define lambda ((coord){0.,0.,0.})
#endif
#endif

// Temporary placement for tangential face gradients

#if EMBED
#define face_avg_gradient_t1_x(a,i)				\
  ((a[1,i-1] + a[1,i] - a[-1,i-1] - a[-1,i])/(4.*Delta))
#define face_avg_gradient_t2_x(a,i) \
  ((a[1,0,i-1] + a[1,0,i] - a[-1,0,i-1] - a[-1,0,i])/(4.*Delta))

#define embed_face_avg_gradient_t1_x(a,i)				\
  (fs.y[]>0. && fs.y[]<1. ? fs.y[1]>0. ? ((a[1,i-1] + a[1,i]) - (a[0,i-1] + a[0,i]))/(2.*Delta) : ((a[0,i-1] + a[0,i])- (a[-1,i] + a[-1,i-1]))/(2.*Delta) \ 
    : fs.y[]>=1. ? (a[1,i-1] + a[1,i] - a[-1,i-1] - a[-1,i])/(4.*Delta) \
    : 0.)
    
#define embed_face_avg_gradient_t2_x(a,i)				\
  (fs.z[]>0. && fs.z[]<1. ? fs.z[1]>0. ? ((a[1,0,i-1] + a[1,0,i]) -(a[0,0,i-1] + a[0,0,i]))/(2.*Delta) : ((a[0,0,i-1] + a[0,0,i])- (a[-1,0,i-1] + a[-1,0,i]))/(2.*Delta) \
    : fs.z[]>=1. ? (a[1,0,i-1] + a[1,0,i] - a[-1,0,i-1] - a[-1,0,i])/(4.*Delta) \
    : 0.)

#define face_avg_gradient_t1_y(a,i) \
  ((a[i-1,1] + a[i,1] - a[i-1,-1] - a[i,-1])/(4.*Delta))
#define face_avg_gradient_t2_y(a,i) \
  ((a[0,1,i-1] + a[0,1,i] - a[0,-1,i-1] - a[0,-1,i])/(4.*Delta))

#define embed_face_avg_gradient_t1_y(a,i)				\
  (fs.x[]>0. && fs.x[]<1. ? fs.x[0,1]>0. ? ((a[i-1,1] + a[i,1]) - (a[i-1,0] + a[i,0]))/(2.*Delta) : ((a[i-1,0] + a[i,0]) - (a[i-1,-1] + a[i,-1]))/(2.*Delta) \ 
      : fs.x[]>=1. ? (a[i-1,1] + a[i,1] - a[i-1,-1] - a[i,-1])/(4.*Delta) \
      : 0.)

#define embed_face_avg_gradient_t2_y(a,i)				\
  (fs.z[]>0. && fs.z[]<1. ? fs.z[0,1]>0. ? ((a[0,1,i-1] + a[0,1,i]) - (a[0,0,i-1] + a[0,0,i]))/(2.*Delta) : ((a[0,0,i-1] + a[0,0,i]) - (a[0,-1,i-1] + a[0,-1,i]))/(2.*Delta) \
  : fs.z[]>=1. ? (a[0,1,i-1] + a[0,1,i] - a[0,-1,i-1] - a[0,-1,i])/(4.*Delta) \
  : 0.)

#define face_avg_gradient_t1_z(a,i) \
  ((a[i-1,0,1] + a[i,0,1] - a[i-1,0,-1] - a[i,0,-1])/(4.*Delta))
#define face_avg_gradient_t2_z(a,i) \
  ((a[0,i-1,1] + a[0,i,1] - a[0,i-1,-1] - a[0,i,-1])/(4.*Delta))

#define embed_face_avg_gradient_t1_z(a,i)				\
  (fs.x[]>0. && fs.x[]<1. ? fs.x[0,0,1]>0. ? ((a[i-1,0,1] + a[i,0,1]) - (a[i-1,0,0] + a[i,0,0]))/(2.*Delta) : ((a[i-1,0,0] + a[i,0,0]) - (a[i-1,0,-1] + a[i,0,-1]))/(2.*Delta) \
    : fs.x[]>=1. ? (a[i-1,0,1] + a[i,0,1] - a[i-1,0,-1] - a[i,0,-1])/(4.*Delta) \
    : 0.)

#define embed_face_avg_gradient_t2_z(a,i)				\
  (fs.y[]>0. && fs.y[]<1. ? fs.y[0,0,1]>0. ? ((a[0,i-1,1] + a[0,i,1]) - (a[0,i-1,0] + a[0,i,0]))/(2.*Delta) : ((a[0,i-1,0] + a[0,i,0]) - (a[0,i-1,-1] + a[0,i,-1]))/(2.*Delta) \
    : fs.y[]>=1. ? (a[0,i-1,1] + a[0,i,1] - a[0,i-1,-1] - a[0,i,-1])/(4.*Delta) \
    : 0.)


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
  double (* embed_flux) (Point, scalar, vector, double *) = p->embed_flux;
#endif // EMBED

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif
  
  boundary({u});
  scalar s=u.x;
  foreach_level_or_leaf (l,nowarning) {

    foreach_dimension() {
      double embedflux=0.;
#if EMBED
      if (embed_flux) {
        if ((cs[]>0.)&&(cs[]<1.)) {
          coord n=facet_normal (point,cs,fs);
          double alpha=plane_alpha(cs[],n);
          coord p;
          double area=plane_area_center(n,alpha,&p);
          if (metric_embed_factor) {
            area*=metric_embed_factor(point,p);
          }
          normalize(&n);
          coord dudn=embed_gradient(point,u,p,n);

          double mua=0.,fua=0.;
          foreach_dimension() {
            mua+=mu.x[]+mu.x[1];
            fua+=fm.x[]+fm.x[1];
          }
          embedflux=-mua/fua*area*(2.*dudn.x*sq(n.x)+dudn.x*sq(n.y)+dudn.y*n.x*n.y)/Delta;
        }
      }
#endif // EMBED
      w.x[] = (dt*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
               #if dimension > 1
		   + mu.y[0,1]*(u.x[0,1] +face_avg_gradient_t1_x (u.y, 1)*Delta)
		   - mu.y[]*(- u.x[0,-1] +
			     face_avg_gradient_t1_x (u.y, 0)*Delta)
               #endif
	       #if dimension > 2
		   + mu.z[0,0,1]*(u.x[0,0,1] +
				  face_avg_gradient_t2_x (u.z, 1)*Delta)
		   - mu.z[]*(- u.x[0,0,-1] +
			     face_avg_gradient_t2_x (u.z, 0)*Delta)
               #endif
       ) + (r.x[] - dt*embedflux)*sq(Delta))/
    (sq(Delta)*(rho[]+lambda.x ) + dt*(2.*mu.x[1] + 2.*mu.x[]
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
  double (* embed_flux) (Point, scalar, vector, double *) = p->embed_flux;
  double (* embed_partial_flux) (Point, scalar, vector, double *,bool,bool) = p->embed_partial_flux;
#endif
    
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  boundary({u});
  scalar s=u.x;
  foreach_dimension() {
    //bool xny=(s.i==u.x.i);
    face vector taux[];
    foreach_face(x)
      taux.x[] = 2.*mu.x[]*face_gradient_x (u.x, 0);
    #if dimension > 1
      foreach_face(y)
	taux.y[] = mu.y[]*(face_gradient_y (u.x, 0) + 
			   embed_face_avg_gradient_t1_x (u.y, 0));
    #endif
    #if dimension > 2
      foreach_face(z)
	taux.z[] = mu.z[]*(face_gradient_z (u.x, 0) + 
			   embed_face_avg_gradient_t2_x (u.z, 0));
    #endif
    foreach (reduction(max:maxres)) {
      double a = 0.;

      double embedflux=0.;
#if EMBED
      if (embed_flux) {
        if ((cs[]>0.)&&(cs[]<1.)) {
          coord n=facet_normal (point,cs,fs);
          double alpha=plane_alpha(cs[],n);
          coord p;
          double area=plane_area_center(n,alpha,&p);
          if (metric_embed_factor) {
            area*=metric_embed_factor(point,p);
          }
          normalize(&n);
          coord dudn=embed_gradient(point,u,p,n);

          double mua=0.,fua=0.;
          foreach_dimension() {
            mua+=mu.x[]+mu.x[1];
            fua+=fm.x[]+fm.x[1];
          }
          embedflux=-mua/fua*area*(2.*dudn.x*sq(n.x)+dudn.x*sq(n.y)+dudn.y*n.x*n.y)/Delta;
        }
      }
#endif // EMBED
      
      foreach_dimension()
	a += taux.x[1] - taux.x[];
      res.x[] = r.x[] - (rho[]+lambda.x)*u.x[] + dt*(a/Delta - embedflux);
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
  boundary (resl);
#else
  /* "naive" discretisation (only 1st order on trees) */
  scalar s=u.x;
  foreach (reduction(max:maxres), nowarning) {
       
    foreach_dimension() {

      double embedflux=0.;
#if EMBED
      if (embed_flux) {
        if ((cs[]>0.)&&(cs[]<1.)) {
          coord n=facet_normal (point,cs,fs);
          double alpha=plane_alpha(cs[],n);
          coord p;
          double area=plane_area_center(n,alpha,&p);
          if (metric_embed_factor) {
            area*=metric_embed_factor(point,p);
          }
          normalize(&n);
          coord dudn=embed_gradient(point,u,p,n);

          double mua=0.,fua=0.;
          foreach_dimension() {
            mua+=mu.x[]+mu.x[1];
            fua+=fm.x[]+fm.x[1];
          }
          embedflux=-mua/fua*area*(2.*dudn.x*sq(n.x)+dudn.x*sq(n.y)+dudn.y*n.x*n.y)/Delta;
        }
      }
#endif // EMBED
      res.x[] = r.x[] - (rho[]+lambda.x)*u.x[] +
        dt*(2.*mu.x[1,0]*face_gradient_x (u.x, 1)
	    - 2.*mu.x[]*face_gradient_x (u.x, 0)
        #if dimension > 1
	    + mu.y[0,1]* (face_gradient_y (u.x, 1) +
			    embed_face_avg_gradient_t1_x (u.y, 1))
	    - mu.y[]*(face_gradient_y (u.x, 0) +
		      embed_face_avg_gradient_t1_x (u.y, 0))
	#endif
        #if dimension > 2
	    + mu.z[0,0,1]*(face_gradient_z (u.x, 1) +
			   embed_face_avg_gradient_t2_x (u.z, 1))
	    - mu.z[]*(face_gradient_z (u.x, 0) +
		      embed_face_avg_gradient_t2_x (u.z, 0))
	#endif
      )/Delta - dt*embedflux;
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
mgstats viscosity (vector u, face vector mu, scalar rho, double dt,
		   int nrelax=4, scalar * res=NULL)
{
  vector r[]; 
  foreach()
    foreach_dimension()
      r.x[] = rho[]*u.x[];

  restriction ({mu, rho});
  struct Viscosity p={mu,rho,dt};

#if EMBED
  p.embed_flux = u.x.boundary[embed] != antisymmetry ? embed_flux : NULL;
#endif // EMBED
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p, nrelax, res,
		   minlevel = 1, // fixme: because of root level
                                 // BGHOSTS = 2 bug on trees
		   tolerance = TOLERANCE_MU? TOLERANCE_MU : TOLERANCE);
}

