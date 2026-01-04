/**
# Shrinking model for phase change simulations
This file implements various shrinking models for the
phase change simulations. The shrinking factor 'zeta' is
computed based on the selected policy and is used in the phase change
model to account for the volume change due to phase change.
This is needed as the phase change in porous media is an ill-posed problem:
the solid mass can vary both in terms of porosity change and in terms of
volume change. 'zeta' is used to saturate the missing degree of freedom.
Of course, this implies direct correlation between the volume change
and the choiche of 'zeta'.
Given that no true physical model exists to determine 'zeta' so far, 
we implement various models to test their effect on the solution.
*/

#include "vofToLs.h"
#include "velocity-potential.h"
#include "common-evaporation.h"

double eps0 = 0.5;
double rhoS = 100.;
double rhoG = 1.;
double muG  = 1e-5;

extern scalar omega;

scalar porosity[];
scalar zeta[];

/**
## Zeta policy
Shrinking policy enumeration.
Used to choose which model to use for the shrinking factor 'zeta'
in the phase change model.
*/

enum zeta_types {
  ZETA_SHRINK,   // Pure shrinking regime
  ZETA_SWELLING, // Fixed interface (no shrinking)
  ZETA_CONST,    // Constant split between shrinking and swelling
  ZETA_SMOOTH,   // Smooth transition between states
  ZETA_SHARP,    // Sharp transition between states
  ZETA_LEVELSET, // Level set method, split based on distance to interface
  ZETA_REACTION  // Split based on the local reaction rate
};

enum zeta_types zeta_policy;

/**
We declare a function to set the shrinking factor 'zeta' based on the
selected policy.
*/

void set_zeta (enum zeta_types zeta_policy) {
  switch (zeta_policy) {
  case ZETA_SHRINK: {
    foreach()
      zeta[] = 1.;
    break;
  }

  case ZETA_SWELLING: {
    foreach()
      zeta[] = 0.;
    break;
  }

  case ZETA_CONST: {
    foreach()
      zeta[] = 0.5;
    break;
  }

  case ZETA_SMOOTH: {
#if AXI
    double radius = pow (3.*statsf(f).sum, 1./3.);
#else
    double radius = sqrt (statsf(f).sum/pi)*2.;
#endif
    foreach()
      zeta[] = 1./(1. + exp(32.*radius - 40.*sqrt(sq(x) + sq(y) + sq(z))));
    break;
  }

  case ZETA_SHARP: {
#if AXI
    double radius = pow (3.*statsf(f).sum, 1./3.);
#else
    double radius = sqrt (statsf(f).sum/pi)*2.;
#endif
    foreach()
      zeta[] = (sqrt (sq(x) + sq(y)) > radius*0.8) ? 1. : 0.;
    break;
  }

  case ZETA_LEVELSET: {
    scalar levelset[];
    vof_to_ls (f, levelset, imax = 5);
    double lmin = statsf(levelset).min;
    if (fabs(lmin) > F_ERR)
      foreach() {
        zeta[] = 1. - levelset[]/statsf(levelset).min;
        zeta[] = clamp(zeta[], 0., 1.);
      }

    vof_to_ls (f, levelset, imax = 60);
    break;
  }

  case ZETA_REACTION: {
    scalar o[];
    foreach()
      o[] = omega[]*f[];

    double o_max = statsf(o).max;

    foreach() {
      zeta[] = o_max > F_ERR ? omega[]/o_max : 0.;
      zeta[] = clamp (zeta[], 0., 1.);
    }
    break;
  }

  default: {
    fprintf (stderr, "Unknown Shrinking model\n");
    return;
  }
  }
}

/**
The porosity is appended as tracer to the volume fraction field 'f'
*/

event defaults (i = 0) {
  f.tracers = list_append (f.tracers, porosity);
  set_zeta (zeta_policy);
}

/**
We reset the gas source term before computing it in the phase change event
*/

event reset_sources (i++) {
  reset ({gas_source}, 0.);
}

event chemistry(i++);

/**
## Phase change event
After the chemistry event, we compute the gas source term based on the 
reaction rate 'omega'. The solid phase velocity field 'ubf' is computed
through the solution of a Poisson equation for the velocity potential 'psi'.
*/

event phasechange(i++) {
  /**
  Clean volume fraction and porosity fields from spurious values.
  */
  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    f[] = (f[] > 1. - F_ERR) ? 1. : f[];
    porosity[] = clamp (porosity[], 0., 1.);
    porosity[] = (f[] > F_ERR) ? porosity[] : 0.;
  }

  set_zeta(zeta_policy);

  mgpsf = project_sv (ubf, psi, fm, mgpsf.nrelax);

  foreach() {
    if (f[] > F_ERR) {
      /**
      porosity is in tracer form, already multiplied by f
      */
#ifdef VARPROP
      gas_source[] = -omega[]*(f[] - porosity[])*(1./rhoGv_S[] - 1./rhoSv[])*cm[];
#else
      gas_source[] = -omega[]*(f[] - porosity[])*(1./rhoG - 1./rhoS)*cm[];
#endif
    }
  }
}

/**
## Interface transport events
We temporarily replace the velocity field 'uf' with the solid phase 
velocity field 'ubf' as this is the actual velocity field used to advect
the interface in the presence of phase change.
*/

face vector ufsave[];
event vof (i++) {
  foreach_face() {
    ufsave.x[] = uf.x[];
    uf.x[] = ubf.x[];
  }
}

/**
We restore the original velocity field after the interface and tracer 
advection has been performed.
*/

event tracer_diffusion (i++) {
  foreach_face()
      uf.x[] = ufsave.x[];
}

/**
The velocity 'ubf' can also be used to determine the timestep.
*/

event stability (i++, last) {
  dt = dtnext (timestep (ubf, dtmax));
}

event cleanup (t = end) {
  free(f.tracers); f.tracers = NULL;
}
