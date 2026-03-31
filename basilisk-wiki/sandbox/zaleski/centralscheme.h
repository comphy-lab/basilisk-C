/**
# Explicit central scheme for the 1D inviscid Burgers equation

Solves
$$
\partial_t u + \partial_x \left(\frac{u^2}{2}\right) = 0
$$
using a first-order forward Euler time step with a centered spatial
difference for the flux divergence:
$$
u^{n+1}_i = u^n_i - \frac{\tau}{2\Delta} \left[ F(u^n_{i+1}) - F(u^n_{i-1}) \right]
$$
where $F(u) = u^2/2$.
*/

#include "run.h"
#include "timestep.h"

/**
## Field declaration
*/

scalar u[];

/**
## Flux function

The Burgers flux $F(u) = u^2/2$.
*/

#define F_burgers(u) (sq(u)/2.)

/**
## Default parameters
*/

event defaults (i = 0)
{
  CFL = 0.4;
}

/**
## Stability and timestep

The CFL condition for Burgers is $\tau \leq \Delta / \max|u|$.
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

Forward Euler with centered flux difference.
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
  Update with centered difference of the flux:
  $$
  u^{n+1}_i = u^n_i - \frac{\tau}{2\Delta}\left[F(u^n_{i+1}) - F(u^n_{i-1})\right]
  $$
  */
  foreach()
    u[] = un[] - (dt/(2.*Delta)) * (F_burgers(un[1]) - F_burgers(un[-1]));
}