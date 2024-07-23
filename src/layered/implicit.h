/**
# Time-implicit barotropic integration

This implements a semi-implicit scheme for the evolution of the
free-surface elevation $\eta$ of the [multilayer solver](README).
The scheme can be summarised as
$$
\begin{aligned}
\frac{\eta^{n + 1} - \eta^n}{\Delta t} & 
  = - \sum_k \nabla \cdot [\theta_H (hu)_k^{n + 1} + (1 - \theta_H)  (hu)^n_k] \\
\frac{(hu)^{n + 1}_k - (hu)_k^n}{\Delta t} &
  =  - \Delta tgh^{n + 1 / 2}_k  (\theta_H \nabla \eta^{n + 1} + 
                                  (1 - \theta_H) \nabla \eta^n)
\end{aligned}
$$
where $\theta_H$ is the "implicitness parameter" typically set to $1/2$. 

The resulting Poisson--Helmholtz equation for $\eta^{n+1}$ is solved
using the multigrid Poisson solver. The convergence statistics are
stored in `mgH`. */

#include "poisson.h"

mgstats mgH;
double theta_H = 0.5;

#define IMPLICIT_H 1

/**
The scheme is unconditionally stable for gravity waves, so the gravity
wave CFL is set to $\infty$, if it has not already been set (typically
by the user). */

event defaults0 (i = 0)
{
  if (CFL_H == 1e40)
    CFL_H = HUGE;
  mgH.nrelax = 4;
}

/**
The relaxation and residual functions of the multigrid solver are
derived from the Poisson--Helmholtz equation for $\eta^{n+1}$ derived
from the equations above
$$
\begin{aligned}
  \eta^{n + 1} + \nabla \cdot (\alpha \nabla \eta^{n + 1}) & = \eta^n -
  \Delta t \sum_k \nabla \cdot (hu)_k^{\star}\\
  \alpha & \equiv - g (\theta \Delta t)^2  \sum_k h^{n + 1 / 2}_k\\
  (hu)_k^{\star} & \equiv (hu)_k^n - \Delta tgh^{n + 1 / 2}_k \theta (1 -
  \theta) \nabla \eta^n
\end{aligned}
$$
*/

trace
static void relax_hydro (scalar * ql, scalar * rhsl, int lev, void * data)
{
  scalar eta = ql[0], rhs_eta = rhsl[0];
  face vector alpha = *((vector *)data);
  foreach_level_or_leaf (lev) {
    double d = - cm[]*Delta;
    double n = d*rhs_eta[];
    eta[] = 0.;
    foreach_dimension() {
      n += alpha.x[0]*a_baro (eta, 0) - alpha.x[1]*a_baro (eta, 1);
      diagonalize (eta) {
	d -= alpha.x[0]*a_baro (eta, 0) - alpha.x[1]*a_baro (eta, 1);
      }
    }
    eta[] = n/d;
  }
}

trace
static double residual_hydro (scalar * ql, scalar * rhsl,
			      scalar * resl, void * data)
{
  scalar eta = ql[0], rhs_eta = rhsl[0], res_eta = resl[0];
  face vector alpha = *((vector *)data);
  double maxres = 0.;
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*a_baro (eta, 0);
  
  foreach (reduction(max:maxres)) {
    res_eta[] = rhs_eta[] - eta[];
    foreach_dimension()
      res_eta[] += (g.x[1] - g.x[])/(Delta*cm[]);
    if (fabs(res_eta[]) > maxres)
      maxres = fabs(res_eta[]);
  }

  return maxres;
}

/**
This can be used to optionally store the residual (for debugging). */

scalar res_eta = {-1};

scalar rhs_eta;
face vector alpha_eta;

/**
The semi-implicit update of the layer heights is done in two
steps. The first step is the explicit advection to time $t + (1 -
\theta_H)\Delta t$ of all tracers (including layer heights) i.e. 
$$
\begin{aligned}
h_k^{n + \theta} & = h_k^n - (1 - \theta_H) \Delta t \nabla \cdot (hu)^n_k
\end{aligned}
$$
*/

event half_advection (i++) {
  if (theta_H < 1.)
    advect (tracers, hu, hf, (1. - theta_H)*dt);
}

/**
The r.h.s. and $\alpha$ coefficients of the Poisson--Helmholtz
equation are computed using the flux values at the "half-timestep". */

event acceleration (i++)
{    
  face vector su[];
  alpha_eta = new face vector;
  double C = - sq(theta_H*dt);
  foreach_face() {
    double ax = theta_H*a_baro (eta, 0);
    su.x[] = 0., alpha_eta.x[] = 0.;
    foreach_layer() {
      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;
      hu.x[] = (1. - theta_H)*(hu.x[] + dt*hf.x[]*ax) + theta_H*hf.x[]*uf;
      hu.x[] += dt*(theta_H*ha.x[] - hf.x[]*ax);
      ha.x[] -= hf.x[]*ax;
      su.x[] += hu.x[];
      alpha_eta.x[] += hf.x[];
    }
    alpha_eta.x[] *= C;
  }

  /**
  The r.h.s. is
  $$
  \text{rhs}_\eta = \eta^n - \Delta t\sum_k\nabla\cdot(hu)^\star_k
  $$
  */
  
  rhs_eta = new scalar;
  foreach() {
    rhs_eta[] = eta[];
    foreach_dimension()
      rhs_eta[] -= dt*(su.x[1] - su.x[])/(Delta*cm[]);
  }

  /**
  The fields used by the relaxation function above (and/or by the
  [relaxation function](nh.h#relax_nh) of the non-hydrostatic solver)
  need to be restricted to all levels. */
  
  restriction ({cm, fm, zb, h, hf, alpha_eta});

  /**
  The restriction function for $\eta$, which has been modified by the
  [multilayer solver](hydro.h#defaults0), needs to be replaced by the
  (default) averaging function for the multigrid solver to work
  properly. */
  
#if TREE
  eta.restriction = restriction_average;
#endif
}

/**
In the second (implicit) step, the Poisson--Helmholtz equation for
$\eta^{n+1}$ is solved and the corresponding values for the fluxes
$(hu)^{n+1}$ are obtained by applying the corresponding pressure
gradient term. */

event pressure (i++)
{
  mgH = mg_solve ({eta}, {rhs_eta}, residual_hydro, relax_hydro, &alpha_eta,
		  res = res_eta.i >= 0 ? (scalar *){res_eta} : NULL,
		  nrelax = 4, minlevel = 1,
		  tolerance = TOLERANCE);
  delete ({rhs_eta, alpha_eta});

  /**
  The restriction function for $\eta$ is restored. */
  
#if TREE
  eta.restriction = restriction_eta;
#endif

  /**
  Note that what is stored in `hu` corresponds to
  $\theta_H(hu)^{n+1}$ since this is the flux which will be used in the
  [pressure event](hydro.h#pressure) to perform the "implicit" update of
  the tracers (including layer heights) i.e. 
  $$
  \begin{aligned}
  h_k^{n + 1} & = h_k^{n + \theta} - \Delta t \nabla \cdot \theta_H (hu)^{n+1}_k
  \end{aligned}
  $$
  */
  
  foreach_face() {
    double ax = theta_H*a_baro (eta, 0);
    foreach_layer() {
      ha.x[] += hf.x[]*ax;
      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;
      hu.x[] = theta_H*(hf.x[]*uf + dt*ha.x[]) - dt*ha.x[];
    }
  }
}

/**
## References

~~~bib
@article{vitousek2013stability,
  title={Stability and consistency of nonhydrostatic free-surface models 
         using the semi-implicit $\theta$-method},
  author={Vitousek, Sean and Fringer, Oliver B},
  journal={International Journal for Numerical Methods in Fluids},
  volume={72},
  number={5},
  pages={550--582},
  year={2013},
  publisher={Wiley Online Library}
}
~~~
*/
