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
 
#. Compute the wresidual $\boldsymbol{res} = \boldsymbol{b} - L(\tilde{\boldsymbol{a}})$.
#. If $\left\lVert\boldsymbol{res}\right\rVert < \epsilon$,
 $\tilde{\boldsymbol{a}}$ is good enough, stop.
#. Else, solve $L(d\boldsymbol{a})\simeq\boldsymbol{res}$.
#. Add $d\boldsymbol{a}$ to $\tilde{\boldsymbol{a}}$ and go back to step 1. 
 
This generic strategy is implemented in the standard Poisson
solver. We also define a data structure for the main parameters of the
viscous problem. */
 
#define INCLUDE_CURVATURE_CORRECTION 1
 
/**
## Selective high-pass Laplace-Beltrami filter
 
When `LB_HP_FILTER` is defined as 1 before including this header the
semi-implicit surface-tension operator is replaced by its high-pass
residual:
$$\mathcal{A}^{HP} = \mathcal{A}^{full} - \mathcal{A}^{smooth}$$
where $\mathcal{A}^{smooth}$ evaluates **every component** of the
surface Laplace-Beltrami â€” diagonal, off-diagonal, and
curvature-gradient â€” with a 5-point (2$\Delta$-spaced) stencil in
place of the standard 3-point ($\Delta$-spaced) stencil.
 
The resulting transfer function is
$$\hat{\mathcal{F}}_{HP}(k\Delta)
  = \hat{\mathcal{F}}_{full}(k\Delta) - \hat{\mathcal{F}}_{smooth}(k\Delta)$$
which gives:
 
- $\hat{\mathcal{F}}_{HP}(\pi)   = 1$ â€” full Nyquist damping **preserved**
- $\hat{\mathcal{F}}_{HP}(0)    = 0$ â€” **zero** dissipation of resolved modes
- $k^4$ roll-off at small $k$   â€” biharmonic (hyper-)viscosity character
 
The stability window is unchanged because the Nyquist mode is passed at
full strength.  All three components of the operator are filtered with the
same wider stencil, so the modified relax and residual operators are
mutually consistent regardless of interface orientation. 

*/

/*

FIXME: It appears that using the FILTERED option for jumps in material properties is necessary in many cases. 
This could potentially be improved.

*/
#ifndef LB_HP_FILTER
# define LB_HP_FILTER 0
#endif
#ifndef DISTANCE_EXTRAPOLATE
# define DISTANCE_EXTRAPOLATE 0
#endif
 
#include "poisson.h"
#include "height_distance.h"
//#include "vof2dist.h"
 
struct Viscosity_with_implicit_ST {
  face vector mu;
  scalar rho;
  double dt;
  vector sigma_n_sq_dirac;
  vector sigma_off_diag_dirac;
  face vector sigma_kappa_gradf;
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
 
static void relax_viscosity_st (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity_with_implicit_ST * p = (struct Viscosity_with_implicit_ST *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  (const) vector sigma_n_sq_dirac = p->sigma_n_sq_dirac;
  (const) vector sigma_off_diag_dirac = p->sigma_off_diag_dirac;
  (const) face vector sigma_kappa_gradf = p->sigma_kappa_gradf;
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
 
#if AXI
  vector axi_corr[];
  foreach_level_or_leaf (l) {
    double sigma_dirac = sigma_n_sq_dirac.x[] + sigma_n_sq_dirac.y[];
    axi_corr.x[] = sq(dt)*sigma_dirac*(u.x[0,1] - u.x[0,-1])*Delta
                   /(2.*max(y, Delta/2.))
#if LB_HP_FILTER
                 /* HP AXI: subtract smooth wide-stencil azimuthal gradient */
                 - sq(dt)*sigma_dirac*(ua.x[0,2] - ua.x[0,-2])*Delta
                   /(4.*max(y, Delta/2.))
#endif
                   ;
    axi_corr.y[] = sq(dt)*sigma_dirac*(u.y[0,1] - u.y[0,-1])*Delta
                   /(2.*max(y, Delta/2.))
#if LB_HP_FILTER
                 - sq(dt)*sigma_dirac*(ua.y[0,2] - ua.y[0,-2])*Delta
                   /(4.*max(y, Delta/2.))
#endif
                   ;
  }
#endif
 
  for (int parity = 0; parity < 2; parity++)
    foreach_level_or_leaf (l)
      if (level == 0 || ((point.i + parity) % 2) != (point.j % 2))
#else
#if dimension > 1
  vector ua = u;
#endif
#if AXI
  vector axi_corr[];
  foreach_level_or_leaf (l) {
    double sigma_dirac = sigma_n_sq_dirac.x[] + sigma_n_sq_dirac.y[];
    axi_corr.x[] = sq(dt)*sigma_dirac*(u.x[0,1] - u.x[0,-1])*Delta
                   /(2.*max(y, Delta/2.))
#if LB_HP_FILTER
                 /* HP AXI: subtract smooth wide-stencil azimuthal gradient */
                 - sq(dt)*sigma_dirac*(ua.x[0,2] - ua.x[0,-2])*Delta
                   /(4.*max(y, Delta/2.))
#endif
                   ;
    axi_corr.y[] = sq(dt)*sigma_dirac*(u.y[0,1] - u.y[0,-1])*Delta
                   /(2.*max(y, Delta/2.))
#if LB_HP_FILTER
                 - sq(dt)*sigma_dirac*(ua.y[0,2] - ua.y[0,-2])*Delta
                   /(4.*max(y, Delta/2.))
#endif
                   ;
  }
#endif
 
  foreach_level_or_leaf (l)
#endif
  {
    foreach_dimension()
      w.x[] = (r.x[]*sq(Delta) 
 
#if AXI
        + axi_corr.x[]
#endif
 
#if LB_HP_FILTER
	/* HP diagonal: full âˆ’ smooth (wide-stencil) Laplacian off-centre values */
	+ sq(dt)*sigma_n_sq_dirac.x[]*((u.x[1]+u.x[-1]) - 0.25*(ua.x[2]+ua.x[-2]))
	+ sq(dt)*sigma_n_sq_dirac.y[]*((u.x[0,1]+u.x[0,-1]) - 0.25*(ua.x[0,2]+ua.x[0,-2]))
#else
	+ sq(dt)*sigma_n_sq_dirac.x[]*(u.x[1]+u.x[-1])
	+ sq(dt)*sigma_n_sq_dirac.y[]*(u.x[0,1]+u.x[0,-1])
#endif
	#if dimension == 3
#if LB_HP_FILTER
	+ sq(dt)*sigma_n_sq_dirac.z[]*((u.x[0,0,1]+u.x[0,0,-1]) - 0.25*(ua.x[0,0,2]+ua.x[0,0,-2]))
#else
	+ sq(dt)*sigma_n_sq_dirac.z[]*(u.x[0,0,1]+u.x[0,0,-1])
#endif
	#endif
	#if dimension == 2
	- sq(dt)/2*sigma_off_diag_dirac.x[]*((ua.x[1,1]-ua.x[-1,1])-(ua.x[1,-1]-ua.x[-1,-1]))
#if LB_HP_FILTER
	/* HP off-diagonal smooth: unconditional 12-point stencil.
	   HP/full = (2-cos kx - cos ky)/4 in [0,1]: no anti-damping at any level.
	   minlevel=1 in mg_solve prevents this stencil from reaching level 0. */
	+ sq(dt)*sigma_off_diag_dirac.x[]*(
	    4*(ua.x[1,1]-ua.x[-1,1]-ua.x[1,-1]+ua.x[-1,-1])
	   +(ua.x[1,2]-ua.x[-1,2]-ua.x[1,-2]+ua.x[-1,-2])
	   +(ua.x[2,1]-ua.x[-2,1]-ua.x[2,-1]+ua.x[-2,-1]))/16.
#endif
#if LB_HP_FILTER
	/* HP curvature-gradient 2D consistent with face-vector residual.
	   Face-vector HP: K*(3*(u[]-u[-1]) - (u[1]-u[-2]))/(4*sq(Delta))
	   Main term:       K_x+*(3u[1])/4  — live, same GS coupling as original.
	   Cross-coupling:  K_x+*ua[-1]/4   — frozen (deferred correction). During
	   the black GS pass ua[-1] would be a fresh red cell; keeping it frozen
	   prevents asymmetric intra-sweep feedback that can push rho(M) > 1.
	   Smooth:          -K_x+*ua[2]/4   — frozen (even offset, must be ua). */
        - sq(dt)/2*4/(fm.x[1]+fm.x[]+fm.y[0,1]+fm.y[])
	  *(fm.x[1]*sigma_kappa_gradf.x[1]*(3.*u.x[1]+u.x[-1]-ua.x[2])/4.
	   -fm.x[]*sigma_kappa_gradf.x[]*(3.*u.x[-1]+u.x[1]-ua.x[-2])/4.
	   +fm.y[0,1]*sigma_kappa_gradf.y[0,1]*(3.*u.x[0,1]+u.x[0,-1]-ua.x[0,2])/4.
	   -fm.y[]*sigma_kappa_gradf.y[]*(3.*u.x[0,-1]+u.x[0,1]-ua.x[0,-2])/4.)
#else
        - sq(dt)/2*4/(fm.x[1]+fm.x[]+fm.y[0,1]+fm.y[])*(fm.x[1]*sigma_kappa_gradf.x[1]*u.x[1]-fm.x[]*sigma_kappa_gradf.x[]*u.x[-1]
							+ fm.y[0,1]*sigma_kappa_gradf.y[0,1]*u.x[0,1] - fm.y[]*sigma_kappa_gradf.y[]*u.x[0,-1])
#endif
	#else
 
	- sq(dt)/2*sigma_off_diag_dirac.z[]*(ua.x[1,1,0]-ua.x[-1,1,0]-ua.x[1,-1,0]+ua.x[-1,-1,0])
	- sq(dt)/2*sigma_off_diag_dirac.y[]*(ua.x[1,0,1]-ua.x[-1,0,1]-ua.x[1,0,-1]+ua.x[-1,0,-1])
	- sq(dt)/2*sigma_off_diag_dirac.x[]*(ua.x[0,1,1]-ua.x[0,-1,1]-ua.x[0,1,-1]+ua.x[0,-1,-1])
#if LB_HP_FILTER
	/* HP 3D off-diagonal (xy pair): unconditional 12-point stencil. */
	+ sq(dt)*sigma_off_diag_dirac.z[]*(
	    4*(ua.x[1,1,0]-ua.x[-1,1,0]-ua.x[1,-1,0]+ua.x[-1,-1,0])
	   +(ua.x[1,2,0]-ua.x[-1,2,0]-ua.x[1,-2,0]+ua.x[-1,-2,0])
	   +(ua.x[2,1,0]-ua.x[-2,1,0]-ua.x[2,-1,0]+ua.x[-2,-1,0]))/16.
	+ sq(dt)*sigma_off_diag_dirac.y[]*(
	    4*(ua.x[1,0,1]-ua.x[-1,0,1]-ua.x[1,0,-1]+ua.x[-1,0,-1])
	   +(ua.x[1,0,2]-ua.x[-1,0,2]-ua.x[1,0,-2]+ua.x[-1,0,-2])
	   +(ua.x[2,0,1]-ua.x[-2,0,1]-ua.x[2,0,-1]+ua.x[-2,0,-1]))/16.
	+ sq(dt)*sigma_off_diag_dirac.x[]*(
	    4*(ua.x[0,1,1]-ua.x[0,-1,1]-ua.x[0,1,-1]+ua.x[0,-1,-1])
	   +(ua.x[0,1,2]-ua.x[0,-1,2]-ua.x[0,1,-2]+ua.x[0,-1,-2])
	   +(ua.x[0,2,1]-ua.x[0,-2,1]-ua.x[0,2,-1]+ua.x[0,-2,-1]))/16.
#endif
#if LB_HP_FILTER
	/* HP curvature-gradient 3D: same frozen-cross-coupling principle as 2D.
	   Main terms 3*u[.] live; cross-coupling ua[opposite] and smooth ua[2x] frozen. */
	- sq(dt)/2*6/(fm.x[1]+fm.x[]+fm.y[0,1]+fm.y[]+fm.z[0,0,1]+fm.z[])
	  *(fm.x[1]*sigma_kappa_gradf.x[1]*(3.*u.x[1]+u.x[-1]-ua.x[2])/4.
	   -fm.x[]*sigma_kappa_gradf.x[]*(3.*u.x[-1]+u.x[1]-ua.x[-2])/4.
	   +fm.y[0,1]*sigma_kappa_gradf.y[0,1]*(3.*u.x[0,1]+u.x[0,-1]-ua.x[0,2])/4.
	   -fm.y[]*sigma_kappa_gradf.y[]*(3.*u.x[0,-1]+u.x[0,1]-ua.x[0,-2])/4.
	   +fm.z[0,0,1]*sigma_kappa_gradf.z[0,0,1]*(3.*u.x[0,0,1]+u.x[0,0,-1]-ua.x[0,0,2])/4.
	   -fm.z[]*sigma_kappa_gradf.z[]*(3.*u.x[0,0,-1]+u.x[0,0,1]-ua.x[0,0,-2])/4.)
#else
	- sq(dt)/2*6/(fm.x[1]+fm.x[]+fm.y[0,1]+fm.y[] + fm.z[0,0,1]+fm.z[])*(fm.x[1]*sigma_kappa_gradf.x[1]*u.x[1]-fm.x[]*sigma_kappa_gradf.x[]*u.x[-1]
                                                        + fm.y[0,1]*sigma_kappa_gradf.y[0,1]*u.x[0,1] - fm.y[]*sigma_kappa_gradf.y[]*u.x[0,-1]
							+ fm.z[0,0,1]*sigma_kappa_gradf.z[0,0,1]*u.x[0,0,1] - fm.z[]*sigma_kappa_gradf.z[]*u.x[0,0,-1])
#endif
 
	#endif
	+ dt/rho[]*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
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
      (lambda.x*sq(Delta)
#if LB_HP_FILTER
	/* HP filter: centre coefficient of u_xx^HP = -(2-1/2) = -3/2 */
	+ 1.5*sq(dt)*(sigma_n_sq_dirac.x[] + sigma_n_sq_dirac.y[]
	#if dimension == 3
	+ sigma_n_sq_dirac.z[]
	#endif
	)
#else
	+ 2*sq(dt)*(sigma_n_sq_dirac.x[] + sigma_n_sq_dirac.y[]
	#if dimension == 3
	+ sigma_n_sq_dirac.z[]
	#endif
	)
#endif
#if LB_HP_FILTER
 #if INCLUDE_CURVATURE_CORRECTION
  + (3./4.)*sq(dt)/2*(
      (fm.x[1]*sigma_kappa_gradf.x[1] + fm.x[]*sigma_kappa_gradf.x[])/(fm.x[1] + fm.x[])
    + (fm.y[0,1]*sigma_kappa_gradf.y[0,1] + fm.y[]*sigma_kappa_gradf.y[])/(fm.y[0,1] + fm.y[])
  #if dimension == 3
    + (fm.z[0,0,1]*sigma_kappa_gradf.z[0,0,1] + fm.z[]*sigma_kappa_gradf.z[])/(fm.z[0,0,1] + fm.z[])
  #endif
  )
  #endif
 
#else
 #if INCLUDE_CURVATURE_CORRECTION
  + sq(dt)/2*(
      (fm.x[1]*sigma_kappa_gradf.x[1] + fm.x[]*sigma_kappa_gradf.x[])/(fm.x[1] + fm.x[])
    + (fm.y[0,1]*sigma_kappa_gradf.y[0,1] + fm.y[]*sigma_kappa_gradf.y[])/(fm.y[0,1] + fm.y[])
  #if dimension == 3
    + (fm.z[0,0,1]*sigma_kappa_gradf.z[0,0,1] + fm.z[]*sigma_kappa_gradf.z[])/(fm.z[0,0,1] + fm.z[])
  #endif
  )
  #endif
#endif
	+ dt/rho[]*(2.*mu.x[1] + 2.*mu.x[]
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
 
static double residual_viscosity_st (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Viscosity_with_implicit_ST * p = (struct Viscosity_with_implicit_ST *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  (const) vector sigma_n_sq_dirac = p->sigma_n_sq_dirac;
  (const) vector sigma_off_diag_dirac = p->sigma_off_diag_dirac;
  (const) face vector sigma_kappa_gradf = p->sigma_kappa_gradf;
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

/*if (!Period.x){
   foreach_boundary(right){
      u.x[2,0] = -u.x[-1,0];
      u.y[2,0] = u.y[-1,0];
   }
   foreach_boundary(left){
      u.x[-2,0] = -u.x[1,0];
      u.y[-2,0] = u.y[1,0];
    }
   }
   if (!Period.y){
   foreach_boundary(top){
      u.x[0,2] = u.x[0,-1];
      u.y[0,2] = -u.y[0,-1];
   }
   foreach_boundary(bottom){
      u.x[0,-2] = u.x[0,1];
      u.y[0,-2] = -u.y[0,1];
    }
   }*/

 
#if AXI
  vector axi_res[];
  foreach() {
    double sigma_dirac = sigma_n_sq_dirac.x[] + sigma_n_sq_dirac.y[];
    axi_res.x[] = sq(dt)*sigma_dirac*(u.x[0,1] - u.x[0,-1])
                  /(2.*Delta*max(y, Delta/2.))
#if LB_HP_FILTER
                /* HP AXI: subtract smooth wide-stencil azimuthal gradient */
                - sq(dt)*sigma_dirac*(u.x[0,2] - u.x[0,-2])
                  /(4.*Delta*max(y, Delta/2.))
#endif
                  ;
    axi_res.y[] = sq(dt)*sigma_dirac*(u.y[0,1] - u.y[0,-1])
                  /(2.*Delta*max(y, Delta/2.))
#if LB_HP_FILTER
                - sq(dt)*sigma_dirac*(u.y[0,2] - u.y[0,-2])
                  /(4.*Delta*max(y, Delta/2.))
#endif
                  ;
  }
#endif
 
 
  foreach_dimension() {
    face vector taux[];
    face vector st_diagx[]; //diagonal component
    face vector st_off_diagx1[]; //off-diagonal component 1
    #if dimension == 3
    face vector st_off_diagx2[]; //off-diagonal component 1
    #endif
    face vector st_curvx[]; //curvature dependent term
    
    foreach_face(x){
      taux.x[] = 2.*mu.x[]*(u.x[] - u.x[-1])/Delta;
#if LB_HP_FILTER
      /* --- HP diagonal (x-face) ---
         Smooth face gradient at face i-1/2: (u[i+1]+u[i]-u[i-1]-u[i-2])/(4dx)
         When differenced in foreach() gives (u[+2]-2u+u[-2])/(4dx^2). */
      st_diagx.x[] = (u.x[] - u.x[-1])/Delta
                   - (u.x[1] + u.x[] - u.x[-1] - u.x[-2])/(4.*Delta);
      /* --- HP off-diagonal (x-face, y-gradient) ---
         Smooth: replace dy-1 stencil with dy-2 stencil, halved coefficient. */
      /* HP off-diagonal (x-face, y-gradient): Type C smooth = (1/2)·standard + (1/2)·wide.
         HP = (2(u[-1,1]+u[0,1]-u[-1,-1]-u[0,-1]) - (u[-1,2]+u[0,2]-u[-1,-2]-u[0,-2]))/(16Δ).
         Cell-level HP: [4(u[1,1]-u[-1,1]-u[1,-1]+u[-1,-1])-(u[1,2]-...)-(u[2,1]-...)]/(16Δ²).
         HP/full ratio = (2-cos kxΔ-cos kyΔ)/4 ∈ [0,1]. */
      st_off_diagx1.x[] = (u.x[-1,1]+u.x[0,1]-u.x[-1,-1]-u.x[0,-1])/8./Delta
                        - (u.x[-1,2]+u.x[0,2]-u.x[-1,-2]-u.x[0,-2])/16./Delta;
      #if dimension == 3
      /* HP off-diagonal (x-face, z-gradient) */
      st_off_diagx2.x[] = (u.x[-1,0,1]+u.x[0,0,1]-u.x[-1,0,-1]-u.x[0,0,-1])/8./Delta
                        - (u.x[-1,0,2]+u.x[0,0,2]-u.x[-1,0,-2]-u.x[0,0,-2])/16./Delta;
      #endif
      /* --- HP curvature-gradient (x-face) ---
         Full:   kn*(u[]-u[-1])/dx^2
         Smooth: kn*(u[1]+u[]-u[-1]-u[-2])/(4*dx^2)
         HP = Full - Smooth = kn*(3(u[]-u[-1])-(u[1]-u[-2]))/(4*dx^2) */
      st_curvx.x[] = sigma_kappa_gradf.x[]
                   *(3.*(u.x[] - u.x[-1]) - (u.x[1] - u.x[-2]))/(4.*sq(Delta));
#else
      st_diagx.x[] = (u.x[] - u.x[-1])/Delta;
      st_off_diagx1.x[] = (u.x[-1,1]+u.x[0,1]-u.x[-1,-1]-u.x[0,-1])/4./Delta;
      #if dimension == 3
      st_off_diagx2.x[] = (u.x[-1,0,1]+u.x[0,0,1]-u.x[-1,0,-1]-u.x[0,0,-1])/4./Delta;
      #endif
      st_curvx.x[] = sigma_kappa_gradf.x[]*(u.x[] - u.x[-1])/sq(Delta);
#endif
    }
    #if dimension > 1
      foreach_face(y){
	taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] +
			   (u.y[1,-1] + u.y[1,0])/4. -
			   (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;
#if LB_HP_FILTER
       /* HP diagonal (y-face): smooth gradient (u[0,1]+u[]-u[0,-1]-u[0,-2])/(4dx) */
       st_diagx.y[] = (u.x[] - u.x[0,-1])/Delta
                    - (u.x[0,1] + u.x[] - u.x[0,-1] - u.x[0,-2])/(4.*Delta);
       /* HP off-diagonal (y-face, x-gradient): wide dx-2 stencil */
       /* HP off-diagonal (y-face, x-gradient): Type C smooth, symmetric analog. */
       st_off_diagx1.y[] = (u.x[1,-1]+u.x[1,0]-u.x[-1,-1]-u.x[-1,0])/8./Delta
                         - (u.x[2,-1]+u.x[2,0]-u.x[-2,-1]-u.x[-2,0])/16./Delta;
       #if dimension == 3
       /* HP off-diagonal (y-face, z-gradient): wide dz-2 stencil */
       st_off_diagx2.y[] = (u.x[0,-1,1]+u.x[0,0,1]-u.x[0,-1,-1]-u.x[0,0,-1])/8./Delta
                         - (u.x[0,-1,2]+u.x[0,0,2]-u.x[0,-1,-2]-u.x[0,0,-2])/16./Delta;
       #endif
       /* HP curvature-gradient (y-face) */
       st_curvx.y[] = sigma_kappa_gradf.y[]
                    *(3.*(u.x[] - u.x[0,-1]) - (u.x[0,1] - u.x[0,-2]))/(4.*sq(Delta));
#else
       st_diagx.y[] = (u.x[] - u.x[0,-1])/Delta;
       st_off_diagx1.y[] = (u.x[1,-1]+u.x[1,0]-u.x[-1,-1]-u.x[-1,0])/4./Delta;
       #if dimension == 3
       st_off_diagx2.y[] = (u.x[0,-1,1]+u.x[0,0,1]-u.x[0,-1,-1]-u.x[0,0,-1])/4./Delta;
       #endif
       st_curvx.y[] = sigma_kappa_gradf.y[]*(u.x[] - u.x[0,-1])/sq(Delta);
#endif
      }
    #endif
    #if dimension > 2
      foreach_face(z){
	taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] +
			   (u.z[1,0,-1] + u.z[1,0,0])/4. -
			   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/Delta;
#if LB_HP_FILTER
	/* HP diagonal (z-face): smooth gradient (u[0,0,1]+u[]-u[0,0,-1]-u[0,0,-2])/(4dx) */
	st_diagx.z[] = (u.x[] - u.x[0,0,-1])/Delta
	             - (u.x[0,0,1] + u.x[] - u.x[0,0,-1] - u.x[0,0,-2])/(4.*Delta);
	/* HP off-diagonal (z-face, x-gradient): wide dx-2 stencil */
	st_off_diagx1.z[] = (u.x[1,0,-1]+u.x[1,0,0]-u.x[-1,0,-1]-u.x[-1,0,0])/8./Delta
	                  - (u.x[2,0,-1]+u.x[2,0,0]-u.x[-2,0,-1]-u.x[-2,0,0])/16./Delta;
	#if dimension == 3
	/* HP off-diagonal (z-face, y-gradient): wide dy-2 stencil */
	st_off_diagx2.z[] = (u.x[0,1,-1]+u.x[0,1,0]-u.x[0,-1,-1]-u.x[0,-1,0])/8./Delta
	                  - (u.x[0,2,-1]+u.x[0,2,0]-u.x[0,-2,-1]-u.x[0,-2,0])/16./Delta;
	#endif
	/* HP curvature-gradient (z-face) */
	st_curvx.z[] = sigma_kappa_gradf.z[]
	             *(3.*(u.x[] - u.x[0,0,-1]) - (u.x[0,0,1] - u.x[0,0,-2]))/(4.*sq(Delta));
#else
	st_diagx.z[] = (u.x[] - u.x[0,0,-1])/Delta;
        st_off_diagx1.z[] = (u.x[1,0,-1]+u.x[1,0,0]-u.x[-1,0,-1]-u.x[-1,0,0])/4./Delta;
        #if dimension == 3
	st_off_diagx2.z[] = (u.x[0,1,-1]+u.x[0,1,0]-u.x[0,-1,-1]-u.x[0,-1,0])/4./Delta;
        #endif
	st_curvx.z[] = sigma_kappa_gradf.z[]*(u.x[] - u.x[0,0,-1])/sq(Delta);
#endif
      }
    #endif
    foreach (reduction(max:maxres)) {
      double d = 0.;
      double d2 = 0.;
 
      foreach_dimension(){
	d += taux.x[1] - taux.x[];
	d2 += sigma_n_sq_dirac.x[]*(st_diagx.x[1] - st_diagx.x[]);
      }
      #if dimension == 2
	coord off_diag = (coord){0,0};
	off_diag.x = (st_off_diagx1.x[1,0] - st_off_diagx1.x[] + st_off_diagx1.y[0,1] - st_off_diagx1.y[]);
 
      #else
	coord off_diag = (coord){0,0,0};
	off_diag.z = (st_off_diagx1.x[1,0,0] - st_off_diagx1.x[] + st_off_diagx1.y[0,1,0] - st_off_diagx1.y[]);
	off_diag.y = (st_off_diagx2.x[1,0,0] - st_off_diagx2.x[] + st_off_diagx1.z[0,0,1] - st_off_diagx1.z[]);
	off_diag.x = (st_off_diagx2.y[0,1,0] - st_off_diagx2.y[0,0,0] + st_off_diagx2.z[0,0,1] - st_off_diagx2.z[]);
	
      #endif 
 
      res.x[] = r.x[] - lambda.x*u.x[] + dt/rho[]*d/Delta
		+ sq(dt)*d2/Delta
		#if dimension == 2
		- sq(dt)*sigma_off_diag_dirac.x[]*(off_diag.x)/Delta
		#else
		- sq(dt)*sigma_off_diag_dirac.x[]*(off_diag.x)/Delta
		- sq(dt)*sigma_off_diag_dirac.y[]*(off_diag.y)/Delta
		- sq(dt)*sigma_off_diag_dirac.z[]*(off_diag.z)/Delta
		#endif
		- sq(dt)*(st_curvx.x[1,0]*fm.x[1,0] + st_curvx.x[]*fm.x[])/(fm.x[1,0] + fm.x[])
		- sq(dt)*(st_curvx.y[0,1]*fm.y[0,1] + st_curvx.y[]*fm.y[])/(fm.y[0,1] + fm.y[])
		#if dimension == 3
		- sq(dt)*(st_curvx.z[0,0,1]*fm.z[0,0,1] + st_curvx.z[]*fm.z[])/(fm.z[0,0,1] + fm.z[])
		#endif
 
#if AXI
	+ axi_res.x[]
#endif
	;
 
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
 
#ifdef PRINT_LB
      /**
      Print coordinates and individual terms of the Laplace-Beltrami operator contribution to the
      residual at interfacial cells, for each velocity component. */
      double interfacial_temp = 0;
      foreach_dimension()
	interfacial_temp += sigma_n_sq_dirac.x[];
      if (interfacial_temp > 0.) {
	double _lb_diag = sq(dt)*d2/Delta;
	double _lb_off_diag = 
	  #if dimension == 2
	  - sq(dt)*sigma_off_diag_dirac.x[]*(off_diag.x)/Delta;
	  #else
	  - sq(dt)*sigma_off_diag_dirac.x[]*(off_diag.x)/Delta
	  - sq(dt)*sigma_off_diag_dirac.y[]*(off_diag.y)/Delta
	  - sq(dt)*sigma_off_diag_dirac.z[]*(off_diag.z)/Delta;
	  #endif
	  double _lb_curv = - sq(dt)*(st_curvx.x[1,0]*fm.x[1,0] + st_curvx.x[]*fm.x[])/(fm.x[1,0] + fm.x[])
	  - sq(dt)*(st_curvx.y[0,1]*fm.y[0,1] + st_curvx.y[]*fm.y[])/(fm.y[0,1] + fm.y[])
	  #if dimension == 3
	  - sq(dt)*(st_curvx.z[0,0,1]*fm.z[0,0,1] + st_curvx.z[]*fm.z[])/(fm.z[0,0,1] + fm.z[])
	  #endif
	  ;
	fprintf (ferr, "LB %12e %12e %12e %12e %12e %12e %12e\n", x, y, z, u.x[], _lb_diag/sq(dt), _lb_off_diag/sq(dt), _lb_curv/sq(dt)*Delta);
      }
#endif
 
 
 
 
 
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
 
//scalar kappa[];
//scalar kappa_test[];
double TOLERANCE_MU = 0.;
 
trace
mgstats viscosity_st (vector u, face vector mu, scalar rho, scalar f, double sigma, double dt,
		   int nrelax = 4, scalar * res = NULL)
{
  
  /**
  The velocity field $\boldsymbol{u}_n$ is provided as an initial
  guess $\tilde{\boldsymbol{a}}$. */
 
  assert (dimension >= 2);
 
  //vertex scalar dist[];
  scalar dist[];
  #if DISTANCE_EXTRAPOLATE
  fprintf(ferr, "extrapolate\n");
  dist[left]   = 3.*dist[]    - 3.*dist[1]      + dist[2];
  dist[right]  = 3.*dist[]    - 3.*dist[-1]     + dist[-2];
  dist[bottom] = 3.*dist[]    - 3.*dist[0, 1]   + dist[0, 2];
  dist[top]    = 3.*dist[]    - 3.*dist[0,-1]   + dist[0,-2];
  # if dimension > 2
  dist[front]  = 3.*dist[]    - 3.*dist[0,0, 1] + dist[0,0, 2];
  dist[back]   = 3.*dist[]    - 3.*dist[0,0,-1] + dist[0,0,-2];
  #endif
  #endif
  scalar kappa[];
 
  vector sigma_n_sq_dirac[];
  vector sigma_off_diag_dirac[];
  face vector sigma_kappa_gradf[];
 
 
  //compute signed distance in narrow band, followed by sigma_n_sq_dirac and sigma_nx_ny_dirac
  //vof2dist(f, dist);
  height_distance(f, dist);
 
  //compute normal by gradient of distance function, then populate squared values and cross term, lastly multiply by sigma norm grad f
  //vector temp1[];
 
  /*
  foreach(){
   foreach_dimension()
        sigma_n_sq_dirac.x[] = 0.0;
   if ((dist[] != nodata) && (dist[0,1] != nodata) && (dist[1,0] != nodata) && (dist[1,1] != nodata)){
     sigma_n_sq_dirac.x[] = (dist[1,0] -  dist[] + dist[1,1] - dist[0,1])/2/Delta;
     sigma_n_sq_dirac.y[] = (dist[0,1] -  dist[] + dist[1,1] - dist[1,0])/2/Delta;
     double nn = sqrt(sq(sigma_n_sq_dirac.x[]) + sq(sigma_n_sq_dirac.y[]));
     sigma_n_sq_dirac.x[] /= nn;
     sigma_n_sq_dirac.y[] /= nn;
     foreach_dimension()
      temp1.x[] = sigma_n_sq_dirac.x[];
   }
    sigma_nx_ny_dirac[] = sigma_n_sq_dirac.x[]*sigma_n_sq_dirac.y[];
    sigma_n_sq_dirac.x[] = sq(sigma_n_sq_dirac.x[]);
    sigma_n_sq_dirac.y[] = sq(sigma_n_sq_dirac.y[]);
 
    double df_dy = (f[1,1] - f[1,-1] + 2*(f[0,1] - f[0,-1]) + f[-1,1] - f[-1,-1])/4/Delta; //unclear if this should be 4 or 8. Depends on defintion of dirac
    double df_dx = (f[1,1] - f[-1,1] + 2*(f[1,0] - f[-1,0]) + f[1,-1] - f[-1,-1])/4/Delta;
    //double dirac = sqrt(sq(df_dx) + sq(df_dy));
    double dirac = max((1.0 - fabs((dist[1,0] +  dist[] + dist[1,1] + dist[0,1])/8/Delta))/Delta, 0);
  
 
    sigma_nx_ny_dirac[] *= sigma*dirac;
    sigma_n_sq_dirac.x[] *= sigma*dirac;
    sigma_n_sq_dirac.y[] *= sigma*dirac;
 
  }*/
 
 foreach(){
   foreach_dimension(){
        sigma_n_sq_dirac.x[] = 0.0;
	sigma_off_diag_dirac.x[] = 0.0;
	//temp1.x[] = 0;
   }
   if (fabs(dist[]) < 2.*L0/((double)(1 << grid->maxdepth))){
     foreach_dimension(){
     sigma_n_sq_dirac.x[] = (dist[1] - dist[-1])/2/Delta;}
     double nn = 0;
     foreach_dimension(){
	nn += sq(sigma_n_sq_dirac.x[]);
     }
     foreach_dimension(){
        sigma_n_sq_dirac.x[] /= sqrt(nn);}
     //foreach_dimension()
     // temp1.x[] = sigma_n_sq_dirac.x[];
 
    #if dimension == 2
    foreach_dimension()
    	sigma_off_diag_dirac.x[] = sigma_n_sq_dirac.x[]*sigma_n_sq_dirac.y[];
    #else
    foreach_dimension()
        sigma_off_diag_dirac.z[] = sigma_n_sq_dirac.x[]*sigma_n_sq_dirac.y[];
 
    #endif
    foreach_dimension()
    	sigma_n_sq_dirac.x[] = 1.-sq(sigma_n_sq_dirac.x[]);
 
    //double df_dy = (f[1,1] - f[1,-1] + 2*(f[0,1] - f[0,-1]) + f[-1,1] - f[-1,-1])/4/Delta; //unclear if this should be 4 or 8. Depends on defintion of dirac
    //double df_dx = (f[1,1] - f[-1,1] + 2*(f[1,0] - f[-1,0]) + f[1,-1] - f[-1,-1])/4/Delta;
    //double dirac = sqrt(sq(df_dx) + sq(df_dy));
    double dirac = max((1.0 - fabs(dist[]/2/Delta))/2/Delta, 0);
    //double dirac = dist[] < 2*Delta ? 0.25*Delta*(1.+cos(M_PI*dist[]/2./Delta)) : 0;
    //double dirac = 0.5/Delta;
    //double epsilon = dist[]/2/Delta;
 
    #ifdef PRINT_LB 
    foreach_dimension(){
      sigma_n_sq_dirac.x[] *= sigma*(0.5*(rho1+rho2));
      sigma_off_diag_dirac.x[] *= sigma*(0.5*(rho1+rho2));}
   
    #else
    foreach_dimension(){
      sigma_n_sq_dirac.x[] *= sigma*dirac/(0.5*(rho1+rho2));
      sigma_off_diag_dirac.x[] *= sigma*dirac/(0.5*(rho1+rho2));}
    #endif
    }
  }
 
  /*foreach(){
    kappa_test[]=((temp1.x[1,0]-temp1.x[-1,0])/Delta/3+(temp1.x[1,1]-temp1.x[-1,1])/Delta/12+(temp1.x[1,-1]-temp1.x[-1,-1])/Delta/12+
                 (temp1.y[0,1]-temp1.y[0,-1])/Delta/3+(temp1.y[1,1]-temp1.y[1,-1])/Delta/12+(temp1.y[-1,1]-temp1.y[-1,-1])/Delta/12);
  }*/
 
  //compute curvature and sigma_kappa_gradf (sigma_kappa_gradf is normalized by delta such that it is mesh independent)
  curvature_hp (f, kappa);
  foreach_face(){
 
	double kappaf =
	  (kappa[] < nodata && kappa[-1] < nodata) ?
	  (kappa[] + kappa[-1])/2. :
	  kappa[] < nodata ? kappa[] :
	  kappa[-1] < nodata ? kappa[-1] :
	  0.;
 
	sigma_kappa_gradf.x[] = sigma*kappaf*(f[] - f[-1])/(0.5*(rho1+rho2));
 
  }
 
  
  vector r[];
  foreach()
    foreach_dimension()
      r.x[] = u.x[];  
 
  /**
  We need $\mu$ and $\rho$ on all levels of the grid. */
  
  restriction ((scalar *){sigma_n_sq_dirac, sigma_off_diag_dirac, sigma_kappa_gradf, fm});
  restriction ({mu,rho});
  struct Viscosity_with_implicit_ST p = { mu, rho, dt, sigma_n_sq_dirac, sigma_off_diag_dirac, sigma_kappa_gradf};
  #ifdef PRINT_LB
  residual_viscosity_st ((scalar *){u}, (scalar *){u}, (scalar *){r}, &p);
  return (mgstats){0};
  #else
  return mg_solve ((scalar *){u}, (scalar *){r}, residual_viscosity_st, relax_viscosity_st, &p, nrelax, res, tolerance = TOLERANCE_MU ? TOLERANCE_MU : TOLERANCE, minlevel = 1);
  //return (mgstats){0};
  #endif
}
 
(const) face vector mu_temp = zerof;
scalar pv[];
vector g2[];
 
pv[right] = neumann (neumann_pressure(ghost));
pv[left]  = neumann (- neumann_pressure(0));
 
#if AXI
                             // is zero on the axis of symmetry
pv[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
pv[top]    = neumann (neumann_pressure(ghost));
pv[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
pv[front]  = neumann (neumann_pressure(ghost));
pv[back]   = neumann (- neumann_pressure(0));
#  endif
#endif
 
scalar psmooth[];
psmooth[left] = neumann(0);
psmooth[right] = neumann(0);
psmooth[top] = neumann(0);
psmooth[bottom] = neumann(0);
#if dimension > 2
psmooth[front] = neumann(0);
psmooth[back] = neumann(0);
#endif 
 
 
mgstats mgpv = {0}, mgps = {0};
vector h_st[];
event defaults (i = 0){
  mgpv = (mgstats){0};
  mgps = (mgstats){0};
  f.height = h_st;
  pv.nodump = psmooth.nodump = true;
}
 
 
/*The current discretisation of curvature requires a 5 point stencil. This is not properly defined at domain boundaries. As a work-around we set 
the values in the second boundary ghost cells manually. The following currently assumes symmetry. */
 
 
event tracer_advection(i++){
   heights(f, f.height);
   boundary({h_st});
   foreach_dimension() {
    h_st.x.stencil.bc = s_centered;
    }
   #if dimension == 2
   if (!Period.x){ 
   foreach_boundary(right)
      h_st.y[2,0] = h_st.y[-1,0];
   foreach_boundary(left)
      h_st.y[-2,0] = h_st.y[1,0];
   }
   if (!Period.y){
   foreach_boundary(top)
      h_st.x[0,2] = h_st.x[0,-1];
   foreach_boundary(bottom)
      h_st.x[0,-2] = h_st.x[0,1];
   }
   #else //dimension == 3 
   if (!Period.x){
   foreach_boundary(right){
      h_st.y[2,0] = h_st.y[-1,0];
      h_st.z[2,0] = h_st.z[-1,0];} 
   foreach_boundary(left){
      h_st.y[-2,0] = h_st.y[1,0];
      h_st.z[-2,0] = h_st.z[1,0];}
   }
   if (!Period.y){
   foreach_boundary(top){
      h_st.x[0,2] = h_st.x[0,-1];
      h_st.z[0,2] = h_st.z[0,-1];}
   foreach_boundary(bottom){
      h_st.x[0,-2] = h_st.x[0,1];
      h_st.z[0,-2] = h_st.z[0,1];}
   }
   if (!Period.z){
   foreach_boundary(front){
      h_st.x[0,0,2] = h_st.x[0,0,-1];
      h_st.y[0,0,2] = h_st.y[0,0,-1];}
   foreach_boundary(back){
      h_st.x[0,0,-2] = h_st.x[0,0,1];
      h_st.y[0,0,-2] = h_st.y[0,0,1];}
   }
   #endif
    
}
 
event advection_term (i++){
 
  /* compute acceleration + pressure gradient at current timestep and add to g*/
    face vector af = a;
    trash ({af});
    foreach_face()
      af.x[] = 0.;
 
 
  event ("acceleration");
  foreach_face()
    uf.x[] = fm.x[]*(dt*a.x[]);
  mgpv = project (uf, pv, alpha, dt, mgpv.nrelax);
  centered_gradient (pv, g2);
  foreach(){
	p[] += pv[]; //we add the newly calculated pv to p. this generally makes the normal poisson-projection step for p converge signficantly quicker
        foreach_dimension()
                g.x[] += g2.x[];
  }
}
 
 
event viscous_term (i++){
 
  correction (dt);
  mgu = viscosity_st (u, mu, rho, f, M_PI/2.*f.sigma, dt, mgu.nrelax);
  //mgu = viscosity_st (u, mu, rho, f, 0.0, dt, mgu.nrelax);
  correction (-dt);
 
  //swap mu with a contant vector such that the code does not use the default viscosity solver 
  swap (face vector, mu, mu_temp);
}
 
 
event end_timestep (i++){
  //swap back 
  swap (face vector, mu, mu_temp);
 
  /*subtract acceleration + pressure gradient at current timestep from g. at this point g represents the pressure gradient that balances the viscous and convective terms*/ 
 
   foreach(){
	psmooth[] = 0;
	p[] -= pv[]; //after projection, we subtract pv from p. if one wishes to visualize the pressure field at the end of the solution loop they should take the sum.
        foreach_dimension()
                g.x[] -= g2.x[];
   }
 
  foreach_face(){
    uf.x[] = fm.x[]*face_value (u.x, 0);
  }
  mgps = project (uf, psmooth, alpha=fm, dt = dt, nrelax = mgps.nrelax);
 
   //DOES THIS IMPROVE RESULTS???
  /*face vector gf[];
  vector g_temp[];
  foreach_face()
    gf.x[] = - fm.x[]*(psmooth[] - psmooth[-1])/Delta;
 
  foreach()
    foreach_dimension()
      g_temp.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
 
  foreach()
    foreach_dimension()
      u.x[] += dt*g_temp.x[];*/
 
  fprintf(ferr, "mgpv.i %d mgpv.nrelax %d mgps.i %d mgps.nrelax %d\n", mgpv.i, mgpv.nrelax, mgps.i, mgps.nrelax);
 
 
}
