/**
# Boiling

Temperature-gradient-driven phase change of pure liquids. The interface is at
the saturation temperature *TIntVal*, and the vaporization rate follows from
the interface energy balance
$$
  \dot{m} = \dfrac{\lambda_l\nabla T_l\cdot\mathbf{n}_\Gamma
    - \lambda_g\nabla T_g\cdot\mathbf{n}_\Gamma}{\Delta h_{ev}}
$$
where $\mathbf{n}_\Gamma$ is the interface normal pointing towards the liquid,
and $\dot{m}$ is positive when the liquid loses mass: the heat conducted from
the superheated phase towards the interface is used for the phase change. This
module includes [phasechange.h](phasechange.h), computes the vaporization
rate and the interface heat fluxes, and sets the phase change model
parameters (`pcm`) for boiling. */

#include "phasechange.h"
#include "gradients.h"

scalar mEvap[], * mEvapList = {mEvap};

/**
## Phase change

The volume and face fractions of each phase are restricted to the fluid
region with embedded boundaries. */

// fixme: put these operations in a function, such that the used can easily move
// it from the phase change event
event phasechange (i++) {
  phase_tracers_to_scalars (liq, f, tol = F_ERR);
  phase_tracers_to_scalars (gas, f, tol = F_ERR);

  scalar fl[], fg[];
  foreach() {
    fl[] = f[];
    fg[] = 1. - f[];
#if EMBED
    fl[] *= cs[];
    fg[] *= cs[];
#endif
  }

  face vector fsl[], fsg[];
  face_fraction (fl, fsl);
  face_fraction (fg, fsg);

  //// fixme as:
  //// face vector fsl[], fsg[];
  //// face_fraction (f, fsl, inverse = false);
  //// face_fraction (f, fsg, inverse = true);

  /**
  The function *plic_gradient()* returns the derivative along the normal
  entering each phase, i.e. $\nabla T_l\cdot\mathbf{n}_\Gamma$ for the liquid
  and $-\nabla T_g\cdot\mathbf{n}_\Gamma$ for the gas: the two contributions
  are therefore summed. The rate is set to zero where the enthalpy of
  evaporation is not positive. */

  foreach_interfacial (f, F_ERR) {
    scalar TL = liq->T, TG = gas->T;
    bool vof = vof_gradients;
    double ltrgrad = plic_gradient (point, TL, fl, fsl, TIntVal, vof, NULL);
    double gtrgrad = plic_gradient (point, TG, fg, fsg, TIntVal, vof, NULL);

    scalar lambdal = liq->lambda, lambdag = gas->lambda;
    scalar deltahev = liq->dhev;

    mEvap[] += (deltahev[] > 0.) ? lambdal[]*ltrgrad/deltahev[] : 0.;
    mEvap[] += (deltahev[] > 0.) ? lambdag[]*gtrgrad/deltahev[] : 0.;
  }

  /**
  The interface conductive heat fluxes enter the temperature equation of each
  phase as explicit source terms. */

  foreach_interfacial_plic (f, F_ERR) {
    scalar TL = liq->T, TG = gas->T;
    bool vof = vof_gradients;
    double ltrgrad = plic_gradient (point, TL, fl, fsl, TIntVal, vof, NULL);
    double gtrgrad = plic_gradient (point, TG, fg, fsg, TIntVal, vof, NULL);

    scalar slT = liq->STexp, sgT = gas->STexp;
    scalar lambdal = liq->lambda, lambdag = gas->lambda;

    slT[] += lambdal[]*ltrgrad*dirac;
    sgT[] += lambdag[]*gtrgrad*dirac;
  }

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);
}

/**
## Phase change model parameters

We set the `pcm` parameters for boiling: the expansion term is shifted to the
gas, `pcm.boiling` selects the boiling treatment of the expansion term, and
the interface regression velocity uses the liquid density. */

// fixme: check compatibility with restore
event init (i = 0) {
  pcm.shifting = SHIFT_TO_GAS;
  pcm.boiling = true;
  pcm.byrhogas = false;
}
