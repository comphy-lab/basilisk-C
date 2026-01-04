/**
# Thermal effects extension for the compressible two-phase solver

This extension of the [two-phase compressible solver](two-phase.h) is
described in detail in [Saade et al, 2023](#saade2023).

This solves the two-phase compressible Navier-Stokes equations
including the total energy equation.
$$
\frac{\partial (f \rho_i)}{\partial t } + 
\nabla \cdot (f \rho_i \mathbf{u}) = 0 
$$
$$
\frac{\partial (\rho_i \mathbf{u})}{\partial t } + 
\nabla \cdot ( \rho_i  \mathbf{u} \mathbf{u}) = 
-\nabla p + \nabla \cdot \tau_i' 
$$
$$
\frac{\partial [\rho_i (e_i + \mathbf{u}^2/2)]}{\partial t } 
+ \nabla \cdot [ \rho_i \mathbf{u} (e_i  + \mathbf{u}^2/2)] =
-\nabla \cdot (\mathbf{u} p_i) + 
\nabla \cdot \left( \tau'_i \mathbf{u} \right)
{\color{blue} - \nabla \cdot \mathbf{q}_i}
$$
an advection equation for the volume fraction $f$
$$
\frac{\partial f}{\partial t} + \mathbf{u} \cdot \nabla f = 0
$$

The term in blue in the energy equation is the heat flux which is
added by this extension. The other terms are treated by the [two-phase
compressible solver](two-phase.h). The heat flux is given by Fourier's
relation as
$$
\nabla \cdot \mathbf{q}_i = - \kappa_i \nabla T_i
$$
where $\kappa_i$ is the thermal conductivity and $T$ is the
temperature field.

The temperature and pressure fields are then coupled through the two
equations (8 and 10 in [Saade et al, 2023](#saade2023))
$$
\rho_ic_{p,i}\frac{DT_i}{Dt} = \beta_iT_i\frac{Dp_i}{Dt} - \nabla\cdot\mathbf{q}_i
$$
with $c_p$ the specific heat capacity and $\beta$ the thermal expansion coefficient.
$$
\left(\frac{\gamma_i}{\rho_ic_i^2} - \frac{\beta_i^2T_i}{\rho_ic_{p,i}}\right)\frac{Dp_i}{Dt} = 
-\frac{\beta_i}{\rho_ic_{p,i}}\nabla\cdot\mathbf{q}_i - \nabla\cdot\mathbf{u}_i
$$
with $\gamma$ the ratio of specific heats.

The system is closed by specifying an Equation Of State (EOS) relating
$p$, $\rho$ and $T$. 

In addition to the primary variables of the compressible two-phase
solver we have the temperature field... */

scalar T[];

/**
... and the thermal conductivities and specific heat capacities. */

double kappa1 = 0., kappa2 = 0.;
double cp1 = 0., cp2 = 0.;

/**
These functions are provided by the Equation Of State. See for example
[the Noble-Able Stiffened-Gas Equation Of State](NASG.h). */

extern double average_temperature (Point point, double p);
extern double thermal_expansion   (Point point);

/**
## Phase-averaged thermal conductivity

By default the arithmetic average of the two phase conductivities is
used to compute the face-centered thermal conductivity. */

const face vector kappa0[] = {0,0,0};
(const) face vector kappa = kappa0;

#ifndef kappa
#define kappa(f) (clamp(f,0.,1.)*(kappa1 - kappa2) + kappa2)
#endif

/**
If the thermal conductivities are not-zero, we need to allocate a new
field. */

event defaults (i = 0)
{
  if (kappa1 || kappa2)
    kappa = new face vector;
}

/**
## Coupled Poisson-Helmholtz equations for temperature and pressure

The robustness of the solver relies on solving equations (8) and (10)
in a strongly-coupled manner with the multigrid solver. After temporal
discretisation this leads to the following coupled Poisson--Helmholtz
system for the temperature $T^{n+1}$ and pressure $p^{n+1}$ (see
eqs. 25 and 26 in [Saade et al, 2023](#saade2023))
$$
\nabla\cdot\left(\kappa\nabla T^{n+1}\right) + \lambda_1 T^{n+1} + \lambda_2 p^{n+1} = \lambda_T
$$
$$
\nabla\cdot\left(\alpha\nabla p^{n+1}\right) + \lambda p^{n+1} + 
\lambda_4\nabla\cdot\left(\kappa\nabla T^{n+1}\right) = \lambda_p
$$
We allocate fields to store the provisionary temperature $T^n$ and the
coefficients. */

scalar Ts[], lambda1[], lambda2[], lambda4[];

/**
We then define the corresponding relaxation and residual functions
which will be used by the [multigrid
solver](/src/poisson.h#mg_solve). */

#include "poisson.h"

struct Thermal {
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};

static void relax_thermal (scalar * al, scalar * bl, int l, void * data)
{
  scalar T = al[0], p = al[1], rhsT = bl[0], rhsp = bl[1];
  struct Thermal * tp = (struct Thermal *) data;
  (const) face vector alpha = tp->alpha;
  (const) scalar lambda = tp->lambda;
  
#if JACOBI
  scalar cT[], cp[];
#else
  scalar cT = T, cp = p;
#endif

  foreach_level_or_leaf (l) {
    double numT = sq(Delta)*(lambda2[]*p[] - rhsT[]), denT = - lambda1[]*sq(Delta);
    double nump = - sq(Delta)*rhsp[], denp = - lambda[]*sq(Delta);
    foreach_dimension() {
      numT += kappa.x[1]*T[1] + kappa.x[]*T[-1];
      denT += kappa.x[1] + kappa.x[];
      nump += alpha.x[1]*p[1] + alpha.x[]*p[-1];
      nump += lambda4[]*(kappa.x[1]*(T[1] - T[]) - kappa.x[]*(T[] - T[-1]));
      denp += alpha.x[1] + alpha.x[];
    }
    cT[] = numT/denT;
    cp[] = nump/denp;
  }

#if JACOBI
  foreach_level_or_leaf (l) {
    T[] = (T[] + 2.*cT[])/3.;
    p[] = (p[] + 2.*cp[])/3.;
  }
#endif
}

static double residual_thermal (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar T = al[0], p = al[1], rhsT = bl[0], rhsp = bl[1], resT = resl[0], resp = resl[1];
  struct Thermal * tp = (struct Thermal *) data;
  (const) face vector alpha = tp->alpha;
  (const) scalar lambda = tp->lambda;
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector gradT[], gradp[];
  foreach_face() {
    gradT.x[] = kappa.x[]*face_gradient_x (T, 0);
    gradp.x[] = alpha.x[]*face_gradient_x (p, 0);
  }
  foreach (reduction(max:maxres)) {
    resT[] = rhsT[] - lambda1[]*T[] - lambda2[]*p[];
    resp[] = rhsp[] - lambda[]*p[];
    foreach_dimension() {
      resT[] -= (gradT.x[1] - gradT.x[])/Delta;
      resp[] -= (gradp.x[1] - gradp.x[])/Delta;
      resp[] -= lambda4[]*(gradT.x[1] - gradT.x[])/Delta;
    }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    resT[] = rhsT[] - lambda1[]*T[] - lambda2[]*p[];
    resp[] = rhsp[] - lambda[]*p[];
    foreach_dimension() {
      resT[] += (kappa.x[0]*face_gradient_x (T, 0) -
		 kappa.x[1]*face_gradient_x (T, 1))/Delta;
      resp[] += (alpha.x[0]*face_gradient_x (p, 0) -
		 alpha.x[1]*face_gradient_x (p, 1))/Delta;
      resp[] += lambda4[]*(kappa.x[0]*face_gradient_x (T, 0) -
			   kappa.x[1]*face_gradient_x (T, 1))/Delta;
    }
#endif // !TREE
    double maxr = max(fabs(resT[]), fabs(resp[]));
    if (maxr > maxres)
      maxres = maxr;
  }
  return maxres;
}

/**
This is the interface for the coupled solver. */
  
mgstats poisson_thermal (scalar p, scalar rhs,
			 (const) face vector alpha,
			 (const) scalar lambda,
			 double tolerance = 0.,
			 int nrelax = 4,
			 int minlevel = 0,
			 scalar * res = NULL)
{
  double defaultol = TOLERANCE;
  if (tolerance)
    TOLERANCE = tolerance;

  restriction ({kappa, alpha, lambda1, lambda2, lambda, lambda4});

  struct Thermal tp = { alpha, lambda, tolerance, nrelax, minlevel, res };
  mgstats s = mg_solve ({T,p}, {Ts,rhs}, residual_thermal, relax_thermal,
			&tp, nrelax, res, minlevel = max(1, minlevel));

  if (tolerance)
    TOLERANCE = defaultol;

  return s;
}

/**
Finally, we overload the *poisson()* function called by [two-phase.h]() with
our new function. */
 
#define poisson(...) poisson_thermal(__VA_ARGS__)
#include "compressible/two-phase.h"
#undef poisson

/**
## Computation of the provisional temperature and Poisson--Helmholtz parameters

This is done before the projection, which will call the coupled Poisson--Helmholtz 
solver. */

event acceleration (i++)
{
  foreach() {
    Ts[] = average_temperature (point, ps[]);
    
    /**
    We compute $\lambda_1$, $\lambda_2$, $\lambda_4$ and r.h.s. for
    the temperature which will be used by the coupled
    Poisson--Helmholtz solver. The $\lambda$ parameter and r.h.s. for
    the pressure are computed by the [two-phase compressible
    solver](two-phase.h#properties). */
    
    double rhocp = cp1*frho1[] + cp2*frho2[];
    lambda1[] = - cm[]*rhocp/dt;
    double beta = thermal_expansion (point);
    lambda2[] = cm[]*beta*Ts[]/dt;
    lambda4[] = beta/rhocp/dt;
    scalar rhsT = Ts;
    rhsT[] = lambda1[]*Ts[] + lambda2[]*ps[];
  }
  
  /**
  We also compute the face-averaged thermal conductivity. */
  
  if (kappa1 || kappa2)
    foreach_face() {
      double ff = (f[] + f[-1])/2.;
      kappa.x[] = fm.x[]*kappa(ff);
    }
}

/**
## Heat flux contribution to the total energy

This adds the $-\nabla\cdot\mathbf{q}_i$ term of the energy
equation. */

event end_timestep (i++)
{
  if (kappa1 || kappa2) {
    face vector gradTv = kappa;
    foreach_face()
      gradTv.x[] = kappa.x[]*face_gradient_x (T,0);
  
    foreach () {
      double energy = 0.;
      foreach_dimension()
	energy += gradTv.x[1] - gradTv.x[];
      energy *= dt/(Delta*cm[]);
      double fc = clamp(f[],0,1);
      fE1[] += energy*fc;
      fE2[] += energy*(1. - fc);
    }
  }
}

/**
## References 

~~~bib
@hal{saade2023, hal-03950917}
~~~
*/
