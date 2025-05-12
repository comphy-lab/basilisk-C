/**
# Temperature Gradient Phase Change Model

This phase change model is suitable for boiling conditions.
The vaporization rate is computed from an energy balance at
the interface (in every interfacial cell), assuming that
the interface is at saturation temperature and that the
heat conduction to the interface is used for the phase change
phenomena. Therefore, the vaporization rate is computed from
the following balance:

$$
\dot{m}_{Evap} \Delta h_{ev} =
    \mathbf{\dot{q}}_l \cdot \mathbf{n_\Gamma}
  + \mathbf{\dot{q}}_g \cdot \mathbf{n_\Gamma} =
  - \lambda_l \left.\dfrac{\partial T_l}{\partial \mathbf{n}_\Gamma}\right\vert_l
  - \lambda_g \left.\dfrac{\partial T_g}{\partial \mathbf{n}_\Gamma}\right\vert_g
$$

where $\lambda$ is the thermal conductivity, $\Delta h_{ev}$
is the enthalpy of evaporation and $\mathbf{n_\Gamma}$ is the
normal to the interface.
After the calculation of the vaporization rate, this model
solves the transport equation for the temperature field $T$
using a two-field formulation and including the presence of
the phase change, both in the diffusion equation and in the
transport of the temperature field which keeps into account
the Stefan convection.
*/

#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"

/**
## Memory Allocations

This phase change model defines a list of vaporization
rates *mEvapList* which is populated by a single scalar
*mEvap* because a single contribution to the vaporization
rate is considered.
*/

#ifndef PHASECHANGE
scalar mEvap[];
scalar * mEvapList = {mEvap};

/**
*fL* and *fG* store the value of volume fractions for the
calculation of the interface gradients, while the vectors
*fsL* and *fsG* contain the face fraction fields computed
using the [fracface.h](/sandbox/ecipriano/src/fracface.h)
module. */

scalar fL[], fG[], f0[];
face vector fsL[], fsG[];

face vector s[];
# define PHASECHANGE
#endif

/**
## Field Allocations

* *T* one-field temperature
* *TL* liquid-phase temperature field
* *TG* gas-phase temperature field
* *TInt* value of temperature at the interface
*/

scalar T[], TL[], TG[], TInt[];

/**
## User Data

Using this phase change model, the user should define the
following variables (SI units):

* *lambda1* Thermal conductivity in liquid phase
* *lambda2* Thermal conductivity in gas phase
* *dhev* Latent heat of vaporization
* *cp1* Specific heat capacity in liquid phase
* *cp2* Specific heat capacity in gas phase
* *TIntVal* Interface temperature value
*/

extern double lambda1, lambda2, dhev, cp1, cp2;
extern double TIntVal;

/**
Allocating diffusivity fields *lambdaf*, volume correction
*theta* and explicit source terms *sgT* and *sTS* for
the diffusion equation. */

face vector lambda1f[], lambda2f[];
scalar thetacorr1[], thetacorr2[];
scalar sgT[], slT[], sgTimp[], slTimp[];

/**
## Defaults

In the defaults event we setup the tracer lists for the
advection of the temperature fields. */

event defaults (i = 0)
{
  TL.inverse = false;
  TG.inverse = true;

#ifdef CONSISTENTPHASE1
  fuext.tracers = list_append (fuext.tracers, TL);
#else
  fu.tracers = list_append (fu.tracers, TL);
#endif
#ifdef CONSISTENTPHASE2
  fuext.tracers = list_append (fuext.tracers, TG);
#else
  fu.tracers = list_append (fu.tracers, TG);
#endif

  /**
  On adaptive meshes, tracers need to use linear interpolation (rather
  than the default bilinear interpolation) to ensure conservation when
  refining cells. */

#if TREE
#if EMBED
      TL.refine = TL.prolongation = refine_embed_linear;
      TG.refine = TG.prolongation = refine_embed_linear;
#else
      TL.refine  = refine_linear;
      TG.refine  = refine_linear;
#endif
      TL.restriction = restriction_volume_average;
      TL.dirty = true; // boundary conditions need to be updated
      TG.restriction = restriction_volume_average;
      TG.dirty = true; // boundary conditions need to be updated
#endif
}

/**
## Init

In the init event, we avoid dumping all the fields that we
don't need to visualize. */

event init (i = 0)
{
  sgT.nodump = true;
  slT.nodump = true;
  fL.nodump  = true;
  fG.nodump  = true;
  thetacorr1.nodump = true;
  thetacorr2.nodump = true;
}

/**
## Finalise

We deallocate the various lists from the memory. */

event cleanup (t = end)
{
  free (fu.tracers), fu.tracers = NULL;
  free (fuext.tracers), fuext.tracers = NULL;
}

/**
## Reset Source Terms

We set to zero phase change source terms, in order to add other
source term contributions from outside this module. */

event reset_sources (i++)
{
  foreach() {
    sgT[] = 0.;
    slT[] = 0.;
    sgTimp[] = 0.;
    slTimp[] = 0.;
  }
}

/**
## Phase Change

In the *phasechange* event, the vaporization rate is computed
and the diffusion step for the mass fraction field (in liquid
and gas phase) is solved. */

event phasechange (i++)
{
  /**
  First we compute the value of the non volume-averaged
  temperature fields. This procedure allows a better
  calculation of the gradients close to the interface. */

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    f0[] = f[];
    fL[] = f[]; fG[] = 1. - f[];
    TL[] = f[] > F_ERR ? TL[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
  }

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  /**
  We compute the vaporization rate from the interface
  jump condition, obtaining the vaporization rate per
  unit of interface surface, stored in *mEvap*. */

  foreach() {
    mEvap[] = 0.; TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      TInt[] = TIntVal;

      double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TIntVal, false);
      double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true,  TIntVal, false);

      mEvap[] += lambda1*ltrgrad/dhev;
      mEvap[] += lambda2*gtrgrad/dhev;

#ifdef SOLVE_LIQONLY
      mEvap[] = lambda1*ltrgrad/dhev;
#endif

#ifdef SOLVE_GASONLY
      mEvap[] = lambda2*gtrgrad/dhev;
#endif
    }
  }

  /**
  The calculation of the interface gradients is used
  also for the calculation of the source terms for the
  diffusion equation of the temperature fields. */

  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fL, fsL), p;
      double alpha = plane_alpha (fL[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TIntVal, false);
      double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true, TIntVal, false);

      double lheatflux = lambda1*ltrgrad;
      double gheatflux = lambda2*gtrgrad;

#ifdef AXI
      slT[] += lheatflux/rho1/cp1*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sgT[] += gheatflux/rho2/cp2*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      slT[] += lheatflux/rho1/cp1*area/Delta*cm[];
      sgT[] += gheatflux/rho2/cp2*area/Delta*cm[];
#endif
    }
  }

  /**
  We restore the tracer form of the liquid and gas-phase
  temperature fields. */

  foreach() {
    TL[] *= f[]*(f[] > F_ERR);
    TG[] *= (1. - f[])*((1. - f[]) > F_ERR);
    T[]  = TL[] + TG[];
  }
}

/**
## Tracer Advection

We let the volume fractions *fu* and *fuext* to
advect the fields YL and YG, as implemented in
the tracer_advection event of [evaporation.h](evaporation.h)
*/

event tracer_advection (i++);

/**
## Tracer Diffusion

We solve the diffusion equations for the temperature fields
accounting for the phase change contributions. */

event tracer_diffusion (i++)
{
  /**
  We remove the fractions of f and mass fractions
  lower than F_ERR and we reconstruct the non-volume
  averaged form of the mass fraction fields, in order
  to improve the discretization of the face gradients
  in the diffusion equation. */

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    //fL[] = f[]; fG[] = 1. - f[];
    //TL[] = fL[] > F_ERR ? TL[]/fL[] : 0.;
    //TG[] = ((1. - fL[]) > F_ERR) ? TG[]/(1. - fL[]) : 0.;

#ifdef CONSISTENTPHASE1
    TL[] = fuext[] > F_ERR ? TL[]/fuext[] : 0.;
#else
    TL[] = fu[] > F_ERR ? TL[]/fu[] : 0.;
#endif
#ifdef CONSISTENTPHASE2
    TG[] = ((1. - fuext[]) > F_ERR) ? TG[]/(1. - fuext[]) : 0.;
#else
    TG[] = ((1. - fu[]) > F_ERR) ? TG[]/(1. - fu[]) : 0.;
#endif
    fL[] = f[]; fG[] = 1. - f[];
  }

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  /**
  We solve the diffusion equations, confined by means of
  the face fraction fields *fsL* and *fsG*. */

  foreach_face() {
    lambda1f.x[] = lambda1/rho1/cp1*fsL.x[]*fm.x[];
    lambda2f.x[] = lambda2/rho2/cp2*fsG.x[]*fm.x[];
  }

  foreach() {
    thetacorr1[] = cm[]*max(fL[], F_ERR);
    thetacorr2[] = cm[]*max(fG[], F_ERR);
    //slT[] = (f[] > F_ERR) ? slT[] : 0.;
    //sgT[] = (f[] > F_ERR) ? sgT[] : 0.;
  }

  diffusion (TG, dt, D=lambda2f, r=sgT, theta=thetacorr2);
  diffusion (TL, dt, D=lambda1f, r=slT, theta=thetacorr1);

  foreach() {
    TL[] *= fL[];
    TG[] *= (1. - fL[]);
  }

  /**
  We reconstruct the one-field temperature field summing
  the two fields $T_L$ and $T_G$ in tracer form. */

  foreach() {
    TL[] = (fL[] > F_ERR) ? TL[] : 0.;
    TG[] = (fG[] > F_ERR) ? TG[] : 0.;
    T[] = TL[] + TG[];
  }
}

