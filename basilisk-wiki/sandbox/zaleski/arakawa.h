/**
# Explicit central skew-symmetric scheme for the 1D inviscid Burgers equation

Solves
$$
\partial_t u + \partial_x \left(\frac{u^2}{2}\right) = 0
$$
using a first-order forward Euler time step with a **skew-symmetric**
centered discretisation of the flux divergence:
$$
u^{n+1}_i = u^n_i - \frac{\Delta t}{6\Delta}
\left[ (u^2_{i+1} - u^2_{i-1}) + u_i(u_{i+1} - u_{i-1}) \right]
$$
This is the average of the conservative form $\partial_x(u^2/2)$ and the
advective form $u\,\partial_x u$, and it preserves the discrete kinetic
$\sum_i u_i^2$ up to $O(\Delta t)^2$ terms before shock formation.
*/

#include "run.h"
#include "timestep.h"

/**
## Field declaration
*/

scalar u[];

/**
## Default parameters
*/

event defaults (i = 0)
{
  CFL = 0.4;
}

/**
## Stability and timestep

The CFL condition for Burgers is $\Delta t \leq \Delta / \max|u|$.
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
## Time integration

Forward Euler with skew-symmetric (energy-conserving) centered difference.

The flux divergence is written as the average of the conservative and
advective forms:
$$
\partial_x\!\left(\frac{u^2}{2}\right)
= \frac{1}{3}\left[\partial_x\!\left(\frac{u^2}{2}\right) + u\,\partial_x u\right]
$$
which discretises as:
$$
u^{n+1}_i = u^n_i
- \frac{\Delta t}{6\Delta}
  \left[(u^2_{i+1} - u^2_{i-1}) + u_i(u_{i+1} - u_{i-1})\right]
$$
One can verify that $\sum_i u_i\,(u^{n+1}_i - u^n_i) = 0$, so
$\sum_i (u^{n+1}_i)^2 = \sum_i (u^n_i)^2$ up to boundary terms.
*/

event time_step (i++, last)
{
  scalar un[];

  /**
  Save the current field before updating.
  */
  foreach()
    un[] = u[];

  /**
  Skew-symmetric update:
  $$
  u^{n+1}_i = u^n_i
  - \frac{\Delta t}{6\Delta}
    \left[(u^2_{i+1} - u^2_{i-1}) + u_i(u_{i+1} - u_{i-1})\right]
  $$
  */
  foreach()
    u[] = un[] - (dt/(6.*Delta)) * (
            (sq(un[1]) - sq(un[-1]))
          + un[]*(un[1] - un[-1])
          );
}