/**
# An implicit solver for the viscous diffusion equation

This header file implicitly solves the viscous diffusion equation:
$$
\rho\frac{\partial\boldsymbol{u}}{\partial t} =
\boldsymbol{\nabla}\cdot\left[\mu\left(\boldsymbol{\nabla u} +
\boldsymbol{\nabla u}^T\right)\right].
$$
Temporally discretised, this equation reads,
$$
-\frac{\rho}{\Delta t}\boldsymbol{u}_{n+1} +
\boldsymbol{\nabla}\cdot\left[\mu\left(\boldsymbol{\nabla u} +
\boldsymbol{\nabla u}^T\right)\right]_{n+1} = -\frac{\rho}{\Delta
t}\boldsymbol{u}_n,
$$
and thus has the form,
$$
L(\boldsymbol{a}) = \boldsymbol{\nabla}\cdot\left[\mu\left(\boldsymbol{\nabla a} +
\boldsymbol{\nabla a}^T\right)\right] = \boldsymbol{b},
$$
where $L()$ is a linear operator, and $\boldsymbol{a}$ and
$\boldsymbol{b}$ are vectors. This system of mutually coupled
equations can therefore be solved efficiently using a multigrid
solver, described for the SGN equations in [Popinet,
2015](/Bibliography#popinet2015). When solving time-dependent
problems, a good initial guess $\tilde{\boldsymbol{a}} =
\boldsymbol{a} - d\boldsymbol{a}$ is available, where
$d\boldsymbol{a}$ is an unknown correction. Therefore, it is usually
more efficient to solve for the equivalent problem,
$$
L(d\boldsymbol{a}) = \boldsymbol{b} - L(\tilde{\boldsymbol{a}}) = \boldsymbol{res},
$$
where $\boldsymbol{res}$ is the residual. Owing to the linearity of
the operator $L()$, $d\boldsymbol{a}$ can be added to the initial
guess $\tilde{\boldsymbol{a}}$, and the process is then repeated until
the residual falls below a given tolerance. The procedure can be
summarised by the following steps:

#. Compute the residual $\boldsymbol{res} = \boldsymbol{b} - L(\tilde{\boldsymbol{a}})$.
#. If $\left\lVert\boldsymbol{res}\right\rVert < \epsilon$,
 $\tilde{\boldsymbol{a}}$ is good enough, stop.
#. Else, solve $L(d\boldsymbol{a})\simeq\boldsymbol{res}$.
#. Add $d\boldsymbol{a}$ to $\tilde{\boldsymbol{a}}$ and go back to step 1. 

This generic strategy is implemented in the standard Poisson
solver. We also define a data structure for the main parameters of the
viscous problem. */

#include "poisson.h"

struct Viscosity {
  face vector mu;
  scalar rho;
  double dt;
};

/**
## Axisymmetry

In Basilisk, axisymmetric simulations are treated in a 2D Cartesian
formulation for code generalisation purposes. The differential element
$dV$ is then accounted for in the metric ([axi.h](/src/axi.h))
incorporated in both $\rho$ and $\mu$. The strain rate tensor
$\boldsymbol{E}$ is of 3D nature:
$$
2\boldsymbol{E} = \boldsymbol{\nabla u} + \boldsymbol{\nabla u}^T =
\begin{bmatrix}
2\frac{\partial u_r}{\partial r} & 0 & \frac{\partial u_r}{\partial z} + 
 \frac{\partial u_z}{\partial r}\\[.1cm]
0 & 2\frac{u_r}{r} & 0\\[.1cm]
\frac{\partial u_r}{\partial z} + \frac{\partial u_z}{\partial r} & 0 & 
 2\frac{\partial u_z}{\partial z}
\end{bmatrix},
$$
so that in cylindrical coordinates, the viscous diffusion equation reads,
$$
\rho\frac{\partial u_r}{\partial t} = \frac{\partial}{\partial r}\left(2\mu\frac{\partial u_r}{\partial r}\right) + \frac{\partial}{\partial z}\left[\mu\left(\frac{\partial u_r}{\partial z} + \frac{\partial u_z}{\partial r}\right)\right] + \frac{2\mu}{r}\left(\frac{\partial u_r}{\partial r} - \frac{u_r}{r}\right),
$$
$$
\rho\frac{\partial u_z}{\partial t} = \frac{\partial}{\partial r}\left[\mu\left(\frac{\partial u_r}{\partial z} + \frac{\partial u_z}{\partial r}\right)\right] + \frac{\partial}{\partial z}\left(2\mu\frac{\partial u_z}{\partial z}\right) + \frac{\mu}{r}\left(\frac{\partial u_r}{\partial z} + \frac{\partial u_z}{\partial r}\right).
$$
Conventionally, in basilisk, $\boldsymbol{e}_y = \boldsymbol{e}_r$ and
$\boldsymbol{e}_x = \boldsymbol{e}_z$. For consistency reasons with 2D
Cartesian simulations, the axisymmetric viscous diffusion equation in
basilisk reads,
$$
\rho'\frac{\partial u}{\partial t} = \frac{\partial}{\partial x}\left(2\mu'\frac{\partial u}{\partial x}\right) + \frac{\partial}{\partial y}\left[\mu'\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)\right],
$$
$$
\rho'\frac{\partial v}{\partial t} = \frac{\partial}{\partial x}\left[\mu'\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)\right] + \frac{\partial}{\partial y}\left(2\mu'\frac{\partial v}{\partial y}\right) - \lambda^*,
$$
where $\rho' = \rho y$ and $\mu' = \mu y$. If one performs the
straightforward derivation, one finds that the equation along $x$ is
identical to the equation along $z$ in cylindrical coordinates. This
is not the case for the equations along $y$ and $r$. That is why the
variable $\lambda^*$ is added to the equation along $y$ in the 2D
Cartesian formulation. Performing the derivations and equating both
equations yields
$$
\lambda^* = \frac{2\mu'}{y^2}v.
$$
The expression under `"lambda.y"` below stems from $\lambda^*$ after
discretisation. In non `AXI` simulations, $\rho' = \rho$, $\mu' = \mu$
and $\lambda^* = 0$. */

#if AXI
# define lambda ((coord){1., 1. + dt/rho[]*(mu.x[] + mu.x[1] + \
					    mu.y[] + mu.y[0,1])/2./sq(y), 0})
#elif SPHERISYM
# define lambda ((coord){1. + 2.*dt/rho[]*(mu.x[] + mu.x[1])/sq(x), 0})
#else // !AXI && !SPHERISYM
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

static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif

  /**
  We have the option of using red/black Gauss-Seidel relaxation or
  "re-use as soon as computed" relaxation. On GPUs (and probably also
  with OpenMP) red/black Gauss-Seidel converges much better (but
  requires two foreach() iterations). Note also that, unlike the other
  option, red/black relaxation should be deterministic. */
  
#if GAUSS_SEIDEL || _GPU
  vector ua[];
  foreach_level (l)
    foreach_dimension()
      ua.x[] = u.x[];
  boundary_level ((scalar *){ua}, l);
  for (int parity = 0; parity < 2; parity++)
    foreach_level_or_leaf (l)
      if (level == 0 || ((point.i + parity) % 2) != (point.j % 2))
#else
#if dimension > 1
  vector ua = u;
#endif
  foreach_level_or_leaf (l)
#endif
  {
    foreach_dimension()
      w.x[] = (r.x[]*sq(Delta) + dt/rho[]*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
#if dimension > 1
					   + mu.y[0,1]*(u.x[0,1] +
							(u.y[1,0] + ua.y[1,1])/4. -
							(u.y[-1,0] + ua.y[-1,1])/4.)
					   - mu.y[]*(- u.x[0,-1] +
						     (ua.y[1,-1] + u.y[1,0])/4. -
						     (ua.y[-1,-1] + u.y[-1,0])/4.)
#endif
#if dimension > 2
					   + mu.z[0,0,1]*(u.x[0,0,1] +
							  (u.z[1,0,0] + ua.z[1,0,1])/4. -
							  (u.z[-1,0,0] + ua.z[-1,0,1])/4.)
					   - mu.z[]*(- u.x[0,0,-1] +
						     (ua.z[1,0,-1] + u.z[1,0,0])/4. -
						     (ua.z[-1,0,-1] + u.z[-1,0,0])/4.)
#endif
					   ))/
      (lambda.x*sq(Delta) + dt/rho[]*(2.*mu.x[1] + 2.*mu.x[]
#if dimension > 1
				      + mu.y[0,1] + mu.y[]
#endif
#if dimension > 2
				      + mu.z[0,0,1] + mu.z[]
#endif
				      ));
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

/**
## Residual computation

This function computes the residual $\boldsymbol{res}$ in step 1 of
the previously described algorithm. It is passed as an argument to the
[multigrid solver](/src/poisson.h#multigrid-solver). */

static double residual_viscosity (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */

  /**
  We manually apply boundary conditions, so that all components are
  treated simultaneously. Otherwise (automatic) BCs would be applied
  component by component before each foreach_face() loop. */
  
  boundary ({u});
  
  foreach_dimension() {
    face vector taux[];
    foreach_face(x)
      taux.x[] = 2.*mu.x[]*(u.x[] - u.x[-1])/Delta;
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
    foreach (reduction(max:maxres)) {
      double d = 0.;
      foreach_dimension()
	d += taux.x[1] - taux.x[];
      res.x[] = r.x[] - lambda.x*u.x[] + dt/rho[]*d/Delta;
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
#else
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres))
    foreach_dimension() {
      res.x[] = r.x[] - lambda.x*u.x[] +
        dt/rho[]*(2.*mu.x[1,0]*(u.x[1] - u.x[])
		  - 2.*mu.x[]*(u.x[] - u.x[-1])
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

#undef lambda

/**
## User interface

A user interface is provided for the solution of the viscous diffusion equation.

### Implicit treatment */

trace
mgstats viscosity (vector u, face vector mu, scalar rho, double dt,
		   int nrelax = 4, scalar * res = NULL)
{
  
  /**
  The velocity field $\boldsymbol{u}_n$ is provided as an initial
  guess $\tilde{\boldsymbol{a}}$. */
  
  vector r[];
  foreach()
    foreach_dimension()
      r.x[] = u.x[];

  /**
  We need $\mu$ and $\rho$ on all levels of the grid. */
  
  restriction ({mu,rho});
  struct Viscosity p = { mu, rho, dt };
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p, nrelax, res);
}

/**
### Explicit treatment

This function does not make use of the multigrid solver. Instead, it
explicitly advances the velocity field in time, in only one
iteration. */

trace
mgstats viscosity_explicit (vector u, face vector mu, scalar rho, double dt)
{
  vector r[];
  mgstats mg = {0};
  struct Viscosity p = { mu, rho, dt };
  mg.resb = residual_viscosity ((scalar *){u}, (scalar *){u}, (scalar *){r}, &p);
  foreach()
    foreach_dimension()
      u.x[] += r.x[];
  return mg;
}

/**
## Possible improvements

#. Note that both functions `residual_viscosity()` and
`relax_viscosity()` are obtained in a very similar manner. The only
difference lies in the formulas of steps 1 and 3, respectively, of the
previously described algorithm. So `relax_viscosity()` is obtained
from `residual_viscosity()` by tweaking and moving the diagonal terms
(which are eventually relaxed) to the left hand side of the equality,
and the residual to its right hand side, or vice versa. It is thus
tedious and unelegant, so a possible improvement is to write a
function encompassing both, where one is automatically derived from
the other, instead of having to copy the respective codes and tweak by
hand. The same goes for the functions written for the
[Poisson-Helmholtz
equation](/src/poisson.h#application-to-the-poissonhelmholtz-equation).
#. Note that these functions are written in "semi-tensorial" form. So a
possible improvement would be to write them in "full tensorial" form,
using nested `foreach_dimension()`. This cannot be done for the moment
as the current version of `foreach_dimension()` does not allow it. */
