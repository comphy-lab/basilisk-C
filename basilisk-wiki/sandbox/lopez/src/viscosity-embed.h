/**
# Viscous solver of the Stokes problem

The Stokes problem appears when the advection term is absent in the
Navier-Stokes equation,
$$
\rho \mathbf{u}_t = \rho \mathbf{g} +
\nabla \cdot [\mu ( \nabla \mathbf{u} + \nabla^T \mathbf{u})]
$$

where $\mathbf{g}$ englobes any superposed acceleration and the
superscript $T$ denotes transverse.

In case of incompressible flow and uniform viscosity, the viscous term
can be further simplified. In effect, since
$$
\nabla \cdot (\mu  \nabla \mathbf{u}) =
(\nabla \mu) \cdot \nabla \mathbf{u} + \mu \nabla (\nabla \cdot  \mathbf{u}) = 0
$$
the viscous term is just,
$$
\nabla \cdot [\mu ( \nabla \mathbf{u} + \nabla^T \mathbf{u})] =
\nabla \cdot (\mu \nabla^T \mathbf{u}) = \nabla^2 (\mu \mathbf{u})
$$
and the velocity components equations decouples each other. For each
component $v_i$, the following discrete implicit equation is solved
using the multigrid solver
$$
\frac{\Delta t} \nabla (\mu \nabla v_i^{n+1}) - (\rho + \lambda_i) v_i^{n+1}
+ \underbrace{\rho (v_i^{n} + g_i\, \Delta t)}_{r_i} = 0
$$
$\lambda_i$ is a possible extra term due to metric. **/

#include "poisson.h"

struct Viscosity {
  vector u;
  face vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
  double (* embed_flux) (Point, scalar, vector, double *);
};

#if AXI
# define lambda ((coord){0, dt*(mu.x[] + mu.x[1]			\
				+ mu.y[] + mu.y[0,1])*sq(cs[])          \
	/(fm.x[] + fm.x[1] + fm.y[] + fm.y[0,1] + SEPS)/(y)})

#else // not AXI
# if dimension == 1
#   define lambda ((coord){0})
# elif dimension == 2
#   define lambda ((coord){0.,0.})
# elif dimension == 3
#   define lambda ((coord){0.,0.,0.})
#endif
#endif

// Note how the relaxation function uses "naive" gradients i.e. not
// the face_gradient_* macros.

static void relax_diffusion (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

  double (* embed_flux) (Point, scalar, vector, double *) = p->embed_flux;
  foreach_level_or_leaf (l) {
    double avgmu = 0.;
    foreach_dimension()
      avgmu += mu.x[] + mu.x[1];
    avgmu = dt*avgmu + SEPS;
    foreach_dimension() {
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.x, mu, &c) : 0.;
      scalar s = u.x;
      double a = 0.;
      foreach_dimension()
	a += mu.x[1]*s[1] + mu.x[]*s[-1];
      u.x[] = (dt*a + (r.x[] - dt*c)*sq(Delta))/
	(sq(Delta)*(rho[] + lambda.x + dt*d) + avgmu);
    }
  }
  
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
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  foreach_dimension() {
    scalar s = u.x;
    face vector g[];
    foreach_face()
      g.x[] = mu.x[]*face_gradient_x (s, 0);
    foreach (reduction(max:maxres), nowarning) {
      double a = 0.;
      foreach_dimension()
	a += g.x[] - g.x[1];
      res.x[] = r.x[] - (rho[] + lambda.x)*u.x[] - dt*a/Delta;
      if (embed_flux) {
	double c, d = embed_flux (point, u.x, mu, &c);
	res.x[] -= dt*(c + d*u.x[]);
      }
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
#else
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres), nowarning)
    foreach_dimension() {
      scalar s = u.x;
      double a = 0.;
      foreach_dimension()
	a += mu.x[0]*face_gradient_x (s, 0) - mu.x[1]*face_gradient_x (s, 1);
      res.x[] = r.x[] - (rho[] + lambda.x)*u.x[] - dt*a/Delta;
      if (embed_flux) {
	double c, d = embed_flux (point, u.x, mu, &c);
	res.x[] -= dt*(c + d*u.x[]);
      }
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
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

  p.embed_flux = u.x.boundary[embed] != antisymmetry ? embed_flux : NULL;
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_diffusion, relax_diffusion, &p, p.nrelax, p.res,
		   minlevel = 1, // fixme: because of root level
                                  // BGHOSTS = 2 bug on trees
		   tolerance = TOLERANCE_MU);
}
