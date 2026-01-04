/**
# Incompressible Navier--Stokes solver (MAC formulation)

We wish to approximate numerically the incompressible Navier--Stokes
equations
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
-\nabla p + \nabla\cdot(\nu\nabla\mathbf{u})
$$
$$
\nabla\cdot\mathbf{u} = 0
$$

We will use the generic time loop, a CFL-limited timestep and we will
need to solve a Poisson problem. */

#include "run.h"
#include "timestep.h"
#include "poisson.h"

/**
The Markers-And-Cells (MAC) formulation was first described in the
pioneering paper of [Harlow and Welch,
1965](/src/references.bib#harlow1965). It relies on a *face*
discretisation of the velocity components `u.x` and `u.y`, relative to
the (centered) pressure `p`. This guarantees the consistency of the
discrete gradient, divergence and Laplacian operators and leads to a
stable (mode-free) integration. */

scalar p[];
face vector u[];

/**
The only parameter is the viscosity coefficient $\nu$.

The statistics for the (multigrid) solution of the Poisson problem are
stored in `mgp`. */

double nu = 0.;
mgstats mgp;

/**
We need to explicitly set the zero normal velocity condition. */

u.n[left]   = 0;
u.n[right]  = 0;
u.n[top]    = 0;
u.n[bottom] = 0;

/**
## Time integration

### Advection--Diffusion 

In a first step, we compute $\mathbf{u}_*$
$$
\frac{\mathbf{u}_* - \mathbf{u}_n}{dt} = \nabla\cdot\mathbf{S}
$$
with $\mathbf{S}$ the symmetric tensor
$$
\mathbf{S} = - \mathbf{u}\otimes\mathbf{u} + \nu\nabla\mathbf{u} =
\left(\begin{array}{cc}
- u_x^2 + 2\nu\partial_xu_x & - u_xu_y + \nu(\partial_yu_x + \partial_xu_y)\\
\ldots & - u_y^2 + 2\nu\partial_yu_y
\end{array}\right)
$$ 

The timestep for this iteration is controlled by the CFL condition
(and the timing of upcoming events). */

double dtmax;

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (u, dtmax));
}

event advance (i++,last)
{

  /**
  We allocate a local symmetric tensor field. To be able to compute the
  divergence of the tensor at the face locations, we need to
  compute the diagonal components at the center of cells and the
  off-diagonal component at the vertices. 
  
  ![Staggering of $\mathbf{u}$ and $\mathbf{S}$](/src/figures/Sxx.svg) */
  
  symmetric tensor S[]; // fixme: this does not work on trees

  /**
  We average the velocity components at the center to compute the
  diagonal components. */

  foreach()
    foreach_dimension()
      S.x.x[] = - sq(u.x[] + u.x[1,0])/4. + 2.*nu*(u.x[1,0] - u.x[])/Delta;

  /**
  We average horizontally and vertically to compute the off-diagonal
  component at the vertices. */

  foreach_vertex()
    S.x.y[] = 
      - (u.x[] + u.x[0,-1])*(u.y[] + u.y[-1,0])/4. +
      nu*(u.x[] - u.x[0,-1] + u.y[] - u.y[-1,0])/Delta;

  /**
  Finally we compute
  $$
  \mathbf{u}_* = \mathbf{u}_n + dt\nabla\cdot\mathbf{S}
  $$ */

  foreach_face()
    u.x[] += dt*(S.x.x[] - S.x.x[-1,0] + S.x.y[0,1] - S.x.y[])/Delta;
}

/**
### Projection 

In a second step we compute
$$
\mathbf{u}_{n+1} = \mathbf{u}_* - \Delta t\nabla p
$$
with the condition
$$
\nabla\cdot\mathbf{u}_{n+1} = 0
$$
This gives the Poisson equation for the pressure
$$
\nabla\cdot(\nabla p) = \frac{\nabla\cdot\mathbf{u}_*}{\Delta t}
$$ */

event projection (i++,last)
{
  const face vector unity[] = {1.[0],1.[0]};
  mgp = project (u, p, unity, dt, mgp.nrelax);
}
