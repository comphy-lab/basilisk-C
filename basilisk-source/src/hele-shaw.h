/**
# Hele-Shaw flow solver

Flows dominated by friction such as *Hele-Shaw flows* or flows in
porous media (*Darcy flows*) can be modelled as
$$
\mathbf{u} = \beta\nabla p
$$
with $p$ analogous to a pressure and where $\beta$ can be a function
of space and time. If the fluid is also incompressible, $p$ needs to
verify the Poisson equation
$$
\nabla\cdot(\beta\nabla p) = \zeta
$$

$\beta$ is often a function of the properties of the fluid such as its
composition, temperature and/or density etc... which also requires the
solution of advection--diffusion equations.

In the following we start from the advection solver and add the
definition of the velocity through the Poisson equation for the
pressure. This is very similar to what is done for the
[streamfunction--vorticity](navier-stokes/stream.h) Navier--Stokes
solver. */

#include "advection.h"
#include "poisson.h"

/**
We allocate the pressure $p$ and divergence field $\zeta$. The
$\beta$ coefficients need to be defined at face locations to
match the locations of the face pressure gradients (and the
face velocity components). These two sets of coefficients are
stored in a vector field. We also allocate space to store the
statistics of the Poisson solver. */

scalar p[], zeta[];
face vector beta[];
mgstats mgp;

/**
We change the default gradient function (used for advection) to
minmod-limited (rather than the centered default). */

event defaults (i = 0)
{
  gradient = minmod2;
}

/**
At every timestep, but after all the other events for this timestep
have been processed (the '`last`' keyword), we update the pressure
field $p$ by solving the Poisson equation with variable coefficient
$\beta$. */

event pressure (i++, last)
{
  mgp = poisson (p, zeta, beta);

  /**
  We then update the velocity field by computing the face pressure
  gradients. */

  trash ({u});
  foreach_face()
    u.x[] = beta.x[]*(p[] - p[-1])/Delta;
}
