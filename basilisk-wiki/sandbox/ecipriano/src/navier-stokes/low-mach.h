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
  \textcolor{blue}{S - \chi_T \dfrac{dP_0}{dt}}
$$
with the deformation tensor $\mathbf{D}=[\nabla\mathbf{u} +
(\nabla\mathbf{u})^T]/2$, $S$ the expansion sources at constant thermodynamic
pressure (temperature, composition, and phase change), and the last term active
only in closed systems (see [low-mach](../../doc/low-mach) for the derivation).
*/

/**
## Field Allocations

We define scalar fields with two possible sources of divergence:
`intexp` is the phase change source term, localized at the
gas-liquid interface, while `drhodt` considers density changes. */

#define LOW_MACH 1

scalar drhodt[], * drhodtlist = NULL;
scalar intexp[], * intexplist = NULL;

bool no_advection_div = false, closed = false;

/**
## Closed systems

If the domain is closed, the total volume cannot change, and the net expansion
must be balanced by a variation of the thermodynamic pressure $P_0$, which is
uniform in space for a low Mach number system. Splitting the density variations
into a thermodynamic pressure contribution and into the remaining effects
(temperature, composition, and phase change), the divergence reads:
$$
  \nabla\cdot\mathbf{u} = S - \chi_T\dfrac{dP_0}{dt}
$$
where $\chi_T = 1/\rho\left(\partial\rho/\partial P\right)_{T,\omega}$ is the
isothermal compressibility. Integrating over the closed domain, where
$\int_\Omega\nabla\cdot\mathbf{u}\,dV = 0$, the pressurization rate is obtained:
$$
  \dfrac{dP_0}{dt} = \dfrac{\int_\Omega S\,dV}{\int_\Omega \chi_T\,dV}
$$
The compensation is weighted on the *local* compressibility. The field `chiT`
must be filled by the module that computes the material properties and if it is
left to zero everywhere, the net expansion is redistributed uniformly over the
domain, which is the behaviour of an incompressible closed system.

If the container exchanges mass through boundaries which are not balanced by an
outlet condition (e.g. an inflow of liquid in a tank whose other boundaries are
walls), the integral of the divergence is equal to the net volumetric flow rate
leaving the domain $\dot{Q} =
\int_{\partial\Omega}\mathbf{u}\cdot\mathbf{n}\,dS$:
$$
  \dfrac{dP_0}{dt} = \dfrac{\int_\Omega S\,dV - \dot{Q}}{\int_\Omega \chi_T\,dV}
$$
The flow rate $\dot{Q}$ is set by the user through the variable `closed_flux`,
which is positive for a net outflow and negative for a net inflow. It is zero by
default, which corresponds to a closed box.

Closed systems are activated with `closed = true`. The user sets the initial
`P0` (null by default) and updates it in the simulation file, e.g.
`P0 += dt*dP0dt`, or from mass conservation, depending on the equation of
state. */

scalar chiT[];
double P0 = 0., dP0dt = 0., _P0 = 0., closed_flux = 0.;

/**
## Virtual metrics

When we use embed, we transport the VOF interface and relative tracers ignoring
the presence of the solid phase (except for the full-solid cells). This
simplifies the VOF transport scheme and avoids small-cells issues. It is
conservative, as long as we ignore the solid interfacial cells when we compute
the conservation. To be consistent with Basilisk, it is useful to transport the
`Phase` quantities in a similar way. Therefore, we define *virtual* metrics,
which allow the classic syntax to be preserved, but adjusting the value of the
metrics. This works also when EMBED and AXI are combined. */

#ifndef cmv
# if EMBED
#   define cmv() (cs[] > 0. ? cm[]/cs[] : 0.)
#   define dvv() (cs[] > 0. ? dv()/cs[] : 0.)
# else
#   define cmv() (cm[])
#   define dvv() (dv())
# endif // EMBED
#endif // cmv

/**
## Projection Function

We define the function that performs the projection step with the volume
expansion term due to the phase change or due to density changes. */

#include "poisson.h"

trace
mgstats project_lowmach (face vector uf, scalar p,
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

  extern int nv;
  static int inv = 0;
  scalar intexp = intexplist[inv];
  scalar drhodt = drhodtlist[inv];

  /**
  The divergence source terms are already multiplied by the metric `cm` (as
  `div[]`), therefore they are integrated over the cell volume without metric,
  while `chiT` is a local property, and it must be multiplied by `cm` to be
  consistent with `div[]`. This makes no difference for Cartesian domains
  without embedded boundaries, where $cm = 1$. */

  double volume = 0., srcsum = 0., chisum = 0.;
  if (closed) {
    foreach (reduction(+:volume) reduction(+:srcsum) reduction(+:chisum)) {
      volume += dvv();
      srcsum += (intexp[] + drhodt[])*dv()/(cm[] + SEPS);
      chisum += chiT[]*dvv();
    }
    /**
    The pressurization rate is computed just once, using the fields of the
    whole domain, because the thermodynamic pressure is a property of the
    system and not of the single velocity field. */

    if (inv == 0)
      dP0dt = (chisum > 0.) ? -(srcsum + closed_flux)/chisum : 0.;
  }

  /**
  If the field `chiT` is zero (by default), the pressurization is evenly
  distributed in the domain. */

  foreach() {
    div[] += (intexp[] + drhodt[])/dt;
    if (closed)
      div[] += cmv()*((chisum > 0.) ? chiT[]*dP0dt : -srcsum/volume)/dt;
  }
  inv++;
  inv = (inv == nv) ? 0 : inv;

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

/**
## Advection Function

We overwrite the advection function for non-incompressible flows, by removing
the divergence of the velocity from the convective term of the momentum
equation, in order to be consistent with the non-conservative formulation.

TODO: embed compatibility */

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
      double fold = f[];
      NOT_UNUSED (fold);
      foreach_dimension() {
        if (no_advection_div)
          f[] += dt*(flux.x[] - flux.x[1] + fold*(u.x[1] - u.x[]))/(Delta*cm[]);
        else
          f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
      }
    }
#else // EMBED
    update_tracer (f, u, flux, dt);
#endif // EMBED
  }

  if (!psrc)
    free (src);
}

/**
We overwrite the functions `project()` and `advection()` in
[centered.h](/src/navier-stokes/centered.h) in order to call `project_lowmach()`
and `advection_div()` instead, accounting for the divergcence source terms. */

#if VELOCITY_JUMP
#define advection(...) advection_div(__VA_ARGS__)
#include "navier-stokes/centered-new.h"
#undef advection
#else
#define project(...) project_lowmach(__VA_ARGS__)
#define advection(...) advection_div(__VA_ARGS__)
#include "navier-stokes/centered-new.h"
#undef advection
#undef project
#endif

/**
We set the default divergence source term to zero (for the liquid phase) */

event defaults (i = 0) {
  dP0dt = 0.;
  _P0 = P0;

  drhodtlist = list_add (drhodtlist, drhodt);
  intexplist = list_add (intexplist, intexp);

  for (int i = 1; i < nv; i++) {
    // Density changes
    {
      scalar drhodt = new scalar;
      char name[80];
      sprintf (name, "drhodt%d", i);
      free (drhodt.name);
      drhodt.name = strdup (name);
      drhodtlist = list_add (drhodtlist, drhodt);
    }

    // Stefan flow
    {
      scalar intexp = new scalar;
      char name[80];
      sprintf (name, "intexp%d", i);
      free (intexp.name);
      intexp.name = strdup (name);
      intexplist = list_add (intexplist, intexp);
    }
  }
}

event cleanup (t = end) {
  P0 = _P0;

  for (int i = 1; i < nv; i++) {
    scalar drhodt = drhodtlist[i];
    scalar intexp = intexplist[i];

    delete ({drhodt,intexp});
  }
  free (drhodtlist), drhodtlist = NULL;
  free (intexplist), intexplist = NULL;
}

