/**
# Incompressible Navier--Stokes solver with divergence source term (centered formulation)

This extension of the [centered.h](/src/navier-stokes/centered.h)
Navier--Stokes equations solver considers situations where the
velocity divergence is not null. Examples are phase change simulations,
or Low Mach systems with variable physical properties.

The modified set of equations reads:

$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} =
  {\color{blue} \dot{m} \left(\dfrac{1}{\rho_g}
  - \dfrac{1}{\rho_l}\right)\delta_\Gamma
  -\color{blue} \dfrac{1}{\rho}\dfrac{D\rho}{Dt}}
$$
with the deformation tensor 
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$.
*/

/**
## Field Allocations

We define scalar fields with two possible sources of divergence:
`stefanflow` is the phase change source term, localized at the
gas-liquid interface, while `drhodt` considers density changes. */

extern scalar stefanflow;
scalar drhodt[];

/**
## Projection Function

We define the function that performs the projection step with the
volume expansion term due to the phase change or due to density
changes. */

#include "poisson.h"

trace
mgstats project_sf (face vector uf, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
  }

  /**
  We add the volume expansion contribution. */

  foreach() {
    div[] += stefanflow[]/dt;
    div[] += drhodt[]/dt;
  }

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson (p, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face()
    uf.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}

#include "utils.h"
#include "bcg.h"

void advection_div (scalar * tracers, face vector u, double dt,
		scalar * src = NULL)
{
  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * psrc = src;
  if (!src)
    for (scalar s in tracers) {
      const scalar zero[] = 0.;
      src = list_append (src, zero);
    }
  assert (list_len (tracers) == list_len (src));

  scalar f, source;
  for (f,source in tracers,src) {
    face vector flux[];
    tracer_fluxes (f, u, flux, dt, source);
#if !EMBED
    foreach() {
#if NO_ADVECTION_DIV
      double fold = f[];
#endif
      foreach_dimension()
#if NO_ADVECTION_DIV
        f[] += dt*(flux.x[] - flux.x[1] + fold*(u.x[1] - u.x[]))/(Delta*cm[]);
#else
        f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
#endif
    }
#else // EMBED
    update_tracer (f, u, flux, dt);
#endif // EMBED
  }

  if (!psrc)
    free (src);
}

/**
We set the default divergence source term to zero (for the liquid phase) */

event defaults (i = 0) {
  foreach()
    drhodt[] = 0.;
}


/**
We overwrite the function `project` in [centered.h](/src/navier-stokes/centered.h)
in order to call `project_sf` instead, accounting for the divergcence
source terms. */

#define project(...) project_sf(__VA_ARGS__)
#define advection(...) advection_div(__VA_ARGS__)
#include "navier-stokes/centered.h"
#undef advection
#undef project
