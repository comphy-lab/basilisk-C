/**
# Azimuthal velocity for axisymmetric flows

The [centered Navier--Stokes solver](centered.h) can be
combined with the [axisymmetric metric](/src/axi.h) but assumes zero
azimuthal velocity ("swirl"). This file adds this azimuthal velocity
component: the $w$ field.

Assuming that $x$ is the axial direction and $y$ the radial direction,
as in [axi.h](/src/axi.h), the incompressible, variable-density and
viscosity, axisymmetric Navier--Stokes equations (with swirl) can be
written
$$
  \partial_x u_x + \partial_y u_y + \frac{u_y}{y} = 0
$$
$$
  \partial_t u_x + u_x \partial_x u_x + u_y \partial_y u_x = -
  \frac{1}{\rho} \partial_x p + \frac{1}{\rho y} \nabla \cdot (2 \mu y \nabla
  \mathbf{D}_x)
$$
$$
  \partial_t u_y + u_x \partial_x u_y + u_y \partial_y u_y - 
  {\color{blue}\frac{w^2}{y}} =
  - \frac{1}{\rho} \partial_y p + \frac{1}{\rho y}  \left( \nabla \cdot (2 \mu
  y \nabla \mathbf{D}_y) - 2 \mu \frac{u_y}{y} \right)
$$
$$
  {\color{blue}
  \partial_t w + u_x \partial_x w + u_y \partial_y w + \frac{u_y w}{y} =
  \frac{1}{\rho y}  \left[ \nabla \cdot (\mu y \nabla w) - w \left(
  \frac{\mu}{y} + \partial_y \mu \right) \right] }
$$
where the terms in blue are the missing "swirl" terms.
We will thus need to solve an advection-diffusion equation for $w$. */

#include "tracer.h"
#include "diffusion.h"

scalar w[], * tracers = {w};

/**
The azimuthal velocity is zero on the axis of symmetry ($y=0$). */

w[bottom] = dirichlet(0);

/**
We will need to add the acceleration term $w^2/y$ in the evolution
equation for $u_y$. If the acceleration field is not allocated yet, we
do so. */

event defaults (i = 0)
{
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
  }
}

/**
The equation for $u_y$ is solved by the centered Navier--Stokes solver
combined with the axisymmetric metric, but the acceleration term
$w^2/y$ is missing. We add it here, taking care of the division by
zero on the axis, and averaging $w$ from cell center to cell face. */

event acceleration (i++)
{
  face vector av = a;
  foreach_face (y)
    av.y[] += y > 0. ? sq(w[] + w[0,-1])/(4.*y) : 0.;
}

/**
The advection of $w$ is done by the [tracer](/src/tracer.h) solver, but we
need to add diffusion. Using the [diffusion](/src/diffusion.h) solver, we
solve
$$
\theta\partial_tw = \nabla\cdot(D\nabla w) + \beta w
$$
Identifying with the diffusion part of the equation for $w$ above, we have
$$
\begin{aligned}
  \theta &= \rho y \\
  D &= \mu y \\
  \beta &= - \left(\rho u_y + \frac{\mu}{y} + \partial_y\mu \right)
\end{aligned}
$$
Note that the *rho* and *mu* fields (defined by the [Navier--Stokes
solver](centered.h)) already include the metric (i.e. are $\rho y$ and
$\mu y$), which explains the divisions by $y$ in the code below. */

event tracer_diffusion (i++)
{
  scalar beta[], theta[];
  foreach() {
    theta[] = rho[];
    double muc = (mu.x[] + mu.x[1] + mu.y[] + mu.y[0,1])/4.;
    double dymu = (mu.y[0,1]/fm.y[0,1] - mu.y[]/fm.y[])/Delta;
    beta[] = - (rho[]*u.y[] + muc/y)/y - dymu;
  }
  diffusion (w, dt, mu, theta = theta, beta = beta);
}
