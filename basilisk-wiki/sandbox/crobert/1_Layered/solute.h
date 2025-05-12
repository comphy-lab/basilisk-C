/** 
# Solutes
The file solute.h aims to add solutes in the multilayer solver.
The solutes will then be advected as any tracers (hydro.h) and diffused 
horizontaly and vertically if the diffusion coefficients are given as attributes in init event.*/
//attribute {
  double D_vert = 0.;
  double D_hor = 0.;
//}

/** By default, a unique solute concentration field $c$ is defined. 
Other solutal fields can be added using tracers and the diffusion coefficient attributes.*/
scalar * solute = NULL;
scalar c;
event defaults (i = 0)
{
  assert (nl > 0);
  c = new scalar[nl];
  reset({c}, 0);

  if (!linearised)
    tracers = list_append (tracers, c);
  solute = list_append (solute, c);
}

/**
# Diffusion
Horizontal and vertical diffusion are called during the viscous event. */
#include "diffusion_solute.h"
event viscous_term (i++)
{
  for (scalar c in solute) {
    if (D_hor > 0)
      horizontal_diffusion (h, c, D_hor, dt);
    boundary ({c});
    if (D_vert > 0)
      foreach()
        vertical_diffusion2 (point, h, c, D_vert, dt);
    boundary ({c});
  }
}


/**
# Evaporation
Evaporation is controlled by the velocity $v_\mathcal{E}$. It must be initialized to add evaporation. 


To ensure stability, the timestep must be maximized. The event set_dtmax is called to ensure that:
$$
dt<\frac{v_e dt}{2}
$$ 
*/
#include "curvature.h"
scalar v_e[];

event stability (i++) 
{
  double dt_max = HUGE;
  foreach()
    if (v_e[] != 0)
      dt_max = min(dt_max, h[0,0,nl-1]/(2*v_e[]));
  dtmax = min(dtmax, dt_max);
}

/**This function calculates the evaporation term.
$$
\frac{dh}{dt} = v_\mathcal{E}\,\mathcal{A}
$$
As the solute does not evaporate, the concentration field must be corrected such as: 
$$
c_{old} * h_{old} = c * h
$$
Evaporation is called just before advection (and eta computation). */

event viscous_term (i++)
{
  foreach()
    if (v_e[] != 0) {
      double h_evap = v_e[] * area_m(point) * dt;
      h[0,0,nl-1] -= h_evap;
      for (scalar c in solute)
        c[0,0,nl-1] += (c[0,0,nl-1] * h_evap/h[0,0,nl-1]);
    }
}