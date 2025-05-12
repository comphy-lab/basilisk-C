/**
# Time-implicit barotropic integration wih surface tension

This implements a semi-implicit scheme for the evolution of the
free-surface elevation $\eta$ of the [multilayer solver](README).
The scheme can be summarised as
$$
\begin{aligned}
\frac{\eta^{n + 1} - \eta^n}{\Delta t} & 
  = - \sum_k \nabla \cdot [\theta_H (hu)_k^{n + 1} + (1 - \theta_H)  (hu)^n_k] \\
\frac{(hu)^{n + 1}_k - (hu)_k^n}{\Delta t} &
  =  - h^{n + 1 / 2}_k  \left(\theta_H \nabla \left(g \eta^{n + 1} - \nabla\cdot\left(\frac{\gamma}{\rho} \frac{\nabla \eta^{n+1}}{\sqrt{1+(\nabla\eta^{n+1})^2}}\right)\right)+ 
  (1 - \theta_H) \nabla \left(g \eta^n - \nabla\cdot\left(\frac{\gamma}{\rho} \frac{\nabla \eta^n}{\sqrt{1+(\nabla\eta^n)^2}}\right)\right)\right)
\end{aligned}
$$
where $\theta_H$ is the "implicitness parameter" typically set to $1/2$. 

The resulting nonlinear Kuramoto--Sivashinsky-like equation for
$\eta^{n+1}$ is solved using the multigrid solver and Newton
iterations. The convergence statistics are stored in `mgH`. */

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
derived from the Kuramoto--Sivashinsky-like equation for $\eta^{n+1}$
derived from the equations above
$$
\begin{aligned}
  \eta^{n + 1} + \nabla \cdot \left(\alpha \nabla \left(g \eta^{n + 1} - \nabla\cdot\left(\frac{\gamma}{\rho}\frac{\nabla \eta^{n+1}}{\sqrt{1+(\nabla\eta^{n+1})^2}}\right)\right)\right) & = \eta^n -
  \Delta t \sum_k \nabla \cdot (hu)_k^{\star}\\
  \alpha & \equiv - (\theta \Delta t)^2  \sum_k h^{n + 1 / 2}_k\\
  (hu)_k^{\star} & \equiv (hu)_k^n - \Delta th^{n + 1 / 2}_k \theta (1 -
  \theta) \nabla \left(g \eta^n - \nabla\cdot\left(\frac{\gamma}{\rho} \frac{\nabla \eta^n}{\sqrt{1+(\nabla\eta^n)^2}}\right)\right)
\end{aligned}
$$
*/

#define cube(x) ((x)*(x)*(x))
#define length(i) (sqrt(1 + sq((eta[] - eta[-1])/Delta)))

struct Newton {
  scalar eta;
  (const) face vector alpha;
};

trace
static void relax_hydrocapillary (scalar * ql, scalar * rhsl, int lev, void * data)
{
  scalar deta = ql[0], rhs_deta = rhsl[0];
  struct Newton * n = (struct Newton *) data;
  (const) face vector alpha = n->alpha;
  scalar eta = n->eta;
  foreach_level_or_leaf (lev) {
    double d = - cm[]*sq(Delta);
    double n = d*rhs_deta[];
    foreach_dimension() {
      n += alpha.x[1]*G*deta[1] + alpha.x[]*G*deta[-1];
      n -= alpha.x[]*sigma[-1]*deta[-2]/(sq(Delta)*cube(length(-1)));
      n += (alpha.x[] * sigma[-1] / (sq(Delta) * cube(length(-1))) + (2.*alpha.x[]+alpha.x[1]) * sigma[] / (sq(Delta) * cube(length(0)))) * deta[-1];
      n += ((alpha.x[]+2.*alpha.x[1]) * sigma[1] / (sq(Delta)*cube(length(1))) + alpha.x[1] * sigma[2] / (sq(Delta)*cube(length(2)))) * deta[1];
      n -= alpha.x[1] * sigma[2]*deta[2]/(sq(Delta)*cube(length(2)));
      d += G*alpha.x[1] + G*alpha.x[];
      d += (2.*alpha.x[]+alpha.x[1]) * sigma[] / (sq(Delta)*cube(length(0)));
      d += (alpha.x[]+2.*alpha.x[1]) * sigma[1] / (sq(Delta)*cube(length(1)));
    }
    deta[] = n/d;
  }
}

trace
static double residual_hydrocapillary (scalar * ql, scalar * rhsl,
			      scalar * resl, void * data)
{
  scalar deta = ql[0], rhs_deta = rhsl[0], res_deta = resl[0];
  struct Newton * n = (struct Newton *) data;
  (const) face vector alpha = n->alpha;
  scalar eta = n->eta;
  double maxres = 0.;
  
  foreach (reduction(max:maxres)) {
    res_deta[] = rhs_deta[] - deta[];
    foreach_dimension() {
      res_deta[] -= -deta[-2] * alpha.x[] * sigma[-1] / (sq(sq(Delta))*cube(length(-1))*cm[]);
      res_deta[] -= deta[-1] * (G*sq(Delta)*alpha.x[]+alpha.x[]*sigma[-1]/cube(length(-1))+(2.*alpha.x[]+alpha.x[1])*sigma[]/cube(length(0)))/(sq(sq(Delta))*cm[]);
      res_deta[] -= deta[] * (-G*alpha.x[]/sq(Delta)-G*alpha.x[1]/sq(Delta)-(2.*alpha.x[]+alpha.x[1])*sigma[]/(sq(sq(Delta))*cube(length(0))) - (alpha.x[]+2.*alpha.x[1])*sigma[1]/(sq(sq(Delta))*cube(length(1)))) / cm[];
      res_deta[] -= deta[1] * (G*sq(Delta)*alpha.x[1]+(alpha.x[]+2.*alpha.x[1])*sigma[1]/cube(length(1))+alpha.x[1]*sigma[2]/cube(length(2)))/(sq(sq(Delta))*cm[]);
      res_deta[] -= -deta[2] * alpha.x[1] * sigma[2] / (sq(sq(Delta))*cube(length(2))*cm[]);
    }
    if (fabs(res_deta[]) > maxres)
      maxres = fabs(res_deta[]);
  }
  boundary (resl);

  return maxres;
}

/**
This can be used to optionally store the residual (for debugging). */

scalar res_deta = {-1};

scalar deta, rhs_eta, rhs_deta;
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
    double ax = - theta_H*gmetric(0)*(G*(eta[] - eta[-1])/Delta + (curvature(0) - curvature(-1))/Delta);
    su.x[] = alpha_eta.x[] = 0.;
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
    alpha_eta.x[] *= C*gmetric(0);
  }
  boundary ((scalar *){alpha_eta, hu, hf, su, ha});

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
  
  // fixme: what about fm?
  restriction ({cm, zb, h, hf, alpha_eta});

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
  // if (t > 1) CFL_H = 50.;
  int newtoniter = 0;
  double newtonerror = HUGE;
  double rhsetanorm, rhsdetanorm;
  struct Newton n;
  n.alpha = alpha_eta;
  n.eta = eta;
  rhs_deta = new scalar;
  deta = new scalar;
  foreach()
    deta[] = 0;
  while (newtonerror > TOLERANCE) {
    newtonerror = 0.;
    newtoniter++;
    foreach() {
      rhs_deta[] = rhs_eta[] - eta[];
      rhs_deta[] -= -alpha_eta.x[]*G*(eta[]-eta[-1])/sq(Delta);
      rhs_deta[] -= alpha_eta.x[]*((eta[-1]-eta[-2])*sigma[-1]/length(-1)-2.*(eta[]-eta[-1])*sigma[]/length(0)+(eta[1]-eta[])*sigma[1]/length(1))/(sq(sq(Delta)));
      rhs_deta[] -= alpha_eta.x[1]*G*(eta[1]-eta[])/sq(Delta);
      rhs_deta[] -= alpha_eta.x[1]*(-(eta[]-eta[-1])*sigma[]/length(0)+2.*(eta[1]-eta[])*sigma[1]/length(1)-(eta[2]-eta[1])*sigma[2]/length(2))/(sq(sq(Delta)));
    }
    boundary({rhs_deta});
    rhsdetanorm = rhsetanorm = 0.;
    foreach() {
      rhsetanorm = max(rhsetanorm, fabs(rhs_eta[]));
      rhsdetanorm = max(rhsdetanorm, fabs(rhs_deta[]));      
    }
    //  fprintf (stderr, "dt = %g ; Current norm of rhs_eta : %g ; current norm of rhs_deta : %g\n", dt, rhsetanorm, rhsdetanorm);
    foreach()
      newtonerror = max(newtonerror, fabs(rhs_deta[]));
    mgH = mg_solve ({deta}, {rhs_deta}, residual_hydrocapillary, relax_hydrocapillary, &n,
                    res = res_deta.i >= 0 ? {res_deta} : NULL,
                    nrelax = 4, minlevel = 1,
                    tolerance = TOLERANCE);
    foreach()
      eta[] += deta[];
    boundary({eta});
    //if (newtoniter == 20) exit(1);
  }
  fprintf (stderr, "t = %g. dt = %g. Newton iteration needed: #%i. Current residual: %g\n", t, dt, newtoniter, newtonerror);
  delete ({deta, rhs_eta, rhs_deta, alpha_eta});

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
    double ax = - theta_H*gmetric(0)*(G*(eta[] - eta[-1])/Delta + (curvature(0) - curvature(-1))/Delta);
    foreach_layer() {
      ha.x[] += hf.x[]*ax;
      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;
      hu.x[] = theta_H*(hf.x[]*uf + dt*ha.x[]) - dt*ha.x[];
    }
  }
  boundary ((scalar *){ha, hu});
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
