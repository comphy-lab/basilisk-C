/**
# A streamfunction--vorticity solver for the Navier--Stokes equations

In two dimensions the incompressible, constant-density Navier--Stokes
equations can be written
$$
\partial_t\omega + \mathbf{u}\cdot\nabla\omega = \nu\nabla^2\omega
$$
$$
\nabla^2\psi = \omega
$$
with $\nu$ the viscosity coefficient. The vorticity $\omega$ and
streamfunction $\psi$ are defined as
$$
\omega = \partial_x u_y - \partial_y u_x
$$
$$
u_x = - \partial_y\psi
$$
$$
u_y = \partial_x\psi
$$
The equation for the vorticity is an advection--diffusion equation
which can be solved using the flux--based advection scheme in
[`advection.h`](/src/advection.h). The equation for the streamfunction
is a [Poisson equation](/src/poisson.h). */

#include "advection.h"
#include "poisson.h"

/**
We allocate the vorticity field $\omega$, the streamfunction field
$\psi$ and a structure to store the statistics on the convergence of
the Poisson solver. The fields advected by the [advection
solver](/src/advection.h) are listed in `tracers`. */

scalar omega[], psi[];
mgstats mgpsi;
scalar * tracers = {omega};

/**
Here we set the default boundary conditions for the
streamfunction. The default convention in Basilisk is no-flow through
the boundaries of the domain, i.e. they are a streamline
i.e. $\psi=$constant on the boundary. */

psi[right]  = dirichlet(0);
psi[left]   = dirichlet(0);
psi[top]    = dirichlet(0);
psi[bottom] = dirichlet(0);

/**
We set the default value for the `CFL` (the default in `utils.h` is
0.5). This is done once at the beginning of the simulation. */

event defaults (i = 0) {
  CFL = 0.8;

  /**
  The default display. */
  
  display ("squares (color = 'omega', spread = -1);");
}

/**
At every timestep we update the streamfunction field $\psi$ by solving
a Poisson equation with the updated vorticity field $\omega$ (which
has just been advected/diffused). */

event velocity (i++)
{
  mgpsi = poisson (psi, omega);

  /**
  Using the new streamfunction, we can then update the components of
  the velocity field. Since they are defined on faces we need to
  average the gradients, which gives the discrete expression below
  (for the horizontal velocity component). The expression for the
  vertical velocity component is obtained by automatic permutation of
  the indices but requires a change of sign: this is done through the
  pseudo-vector `f`. */

  trash ({u});
  coord f = {-1.[0],1.[0]};
  foreach_face()
    u.x[] = f.x*(psi[0,1] + psi[-1,1] - psi[0,-1] - psi[-1,-1])/(4.*Delta);
}
