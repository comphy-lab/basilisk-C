/**
# Skew-symmetric scheme with iterated midpoint time stepping
for the 1D inviscid Burgers equation

Solves
$$
\partial_t u + \partial_x \left(\frac{u^2}{2}\right) = 0
$$
using an **iterated midpoint** time step with a **skew-symmetric**
centered discretisation of the flux divergence.

The skew-symmetric flux divergence operator is:
$$
L(u)_i = \frac{1}{6\Delta}
\left[ (u^2_{i+1} - u^2_{i-1}) + u_i(u_{i+1} - u_{i-1}) \right]
$$

The iterated midpoint scheme reads, starting from $u^{(0)} = u^n$:
$$
u^{(k+1)}_i = u^n_i - \Delta t \, L\!\left(\frac{u^n + u^{(k)}}{2}\right)_i
\quad k = 0, 1, \ldots, N_{\rm iter}-1
$$
and setting $u^{n+1} = u^{(N_{\rm iter})}$.

- $N_{\rm iter} = 1$ recovers the forward Euler scheme
- $N_{\rm iter} = 2$ recovers the RK2 midpoint scheme
- $N_{\rm iter} = 3, 4, \ldots$ converge toward the exactly energy-conserving
  implicit midpoint rule

The number of iterations is controlled by the user-settable parameter
`NITER` (default 3).

The residual of the fixed-point iteration
$$
r^{(k)} = \max_i \left| u^{(k+1)}_i - u^{(k)}_i \right|
$$
is written to the file `residual` with columns `t k r`.
*/

#include "run.h"
#include "timestep.h"

/**
## Field declaration
*/

scalar u[];

/**
## Number of iterations (user can override before including this file)
*/

#ifndef NITER
# define NITER 3
#endif

/**
## Skew-symmetric flux divergence operator
*/

static inline double skew (scalar v, Point point)
{
  return (sq(v[1]) - sq(v[-1]) + v[]*(v[1] - v[-1])) / 6.;
}

/**
## Skew-symmetric flux divergence of the midpoint average $\frac{u^n+u^{(k)}}{2}$
*/

static inline double skew_mid (scalar vn, scalar vk, Point point)
{
  double vm  = 0.5*(vn[]   + vk[]);
  double vmp = 0.5*(vn[1]  + vk[1]);
  double vmm = 0.5*(vn[-1] + vk[-1]);
  return (sq(vmp) - sq(vmm) + vm*(vmp - vmm)) / 6.;
}

/**
## Default parameters
*/

event defaults (i = 0)
{
  CFL = 0.4;
}

/**
## Stability and timestep
*/

double dtmax;

event init (i = 0)
{
  dtmax = DT;
  event ("stability");
}

event set_dtmax (i++, last) dtmax = DT;

event stability (i++, last)
{
  vector vel[];
  foreach()
    vel.x[] = u[];
  dt = dtnext (timestep (vel, dtmax));
}

/**
## Time integration: iterated midpoint scheme
*/

event time_step (i++, last)
{
  scalar un[], uk[], ukp[];
  static FILE * fpr = NULL;

  /**
  Open the residual file on first call.
  */
  if (fpr == NULL)
    fpr = fopen ("residual", "w");

  /**
  Save $u^n$ and initialise the iterate $u^{(0)} = u^n$.
  */
  foreach()
    un[] = uk[] = u[];

  /**
  ### Fixed-point iteration
  $$
  u^{(k+1)}_i = u^n_i - \frac{\Delta t}{\Delta}
  L\!\left(\frac{u^n + u^{(k)}}{2}\right)_i
  $$
  We compute the residual $r^{(k)} = \max_i |u^{(k+1)}_i - u^{(k)}_i|$
  at each iteration and write it to the residual file.
  */
  for (int k = 0; k < NITER; k++) {
    double res = 0.;
    foreach() {
      ukp[] = un[] - (dt/Delta) * skew_mid(un, uk, point);
      res = max(res, fabs(ukp[] - uk[]));
    }
    if(k==NITER-1) {
      foreach()
        fprintf (fpr, "%g %g \n",x,ukp[] - uk[]);
      fprintf (fpr, "\n\n");
    }
    foreach()
      uk[] = ukp[];
  }

  /**
  Set $u^{n+1} = u^{(N_{\rm iter})}$.
  */
  foreach()
    u[] = uk[];
}
