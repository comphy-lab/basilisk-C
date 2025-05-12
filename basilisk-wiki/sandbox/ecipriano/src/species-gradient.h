/**
# Species Gradient Phase Change Model

This phase change model is suitable for vaporization driven
by a gradient of chemical species, between the interface
and the sorrounding gas phase. The vaporization rate (per
unit of interface surface) is computed from a mass balance
at the interface, assuming that the mass fraction on the gas
phase side is constant (i.e. pure liquid phase and constant
vapor pressure). Therefore, the vaporization rate is computed
from the following balance:

$$
\dot{m}_{Evap} = -\dfrac{\rho_g\mathcal{D}_g}{1-\hat{Y}_G}
\left.\dfrac{\partial Y}{\partial \mathbf{n}_\Gamma}\right\vert_g
$$

where $\mathcal{D}_G$ is the diffusivity coefficient in the gas phase,
while $\hat{Y}_G$ is the saturation mass fraction on the
gas phase side. If *DIFFUSIVE* conditions are defined, the
stefan convection is neglected from the mass balance, and
the vaporization rate is computed with the same formula
but neglecting the denominator:

$$
\dot{m}_{Evap} = -\rho_g\mathcal{D}_g
\left.\dfrac{\partial Y}{\partial \mathbf{n}_\Gamma}\right\vert_g
$$

which results in a vaporization rate proportional to the diffusive
flux in the gas phase.
After the calculation of the vaporization rate, this model
solves the transport equation for the mass fraction field $Y$
using a two-field formulation and including the presence of
the phase change, both in the diffusion equation and in the
transport of the mass fraction field which keeps into account
the Stefan convection.
*/

#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"

/**
## Memory Allocations

This phase change model defines a list of vaporization
rates *mEvapList* which is populated by a single scalar
*mEvap* because a single chemical species in liquid phase is
considered.
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

# define PHASECHANGE
#endif

/**
## Field Allocations

* *Y* one-field mass fraction of the chemical species to be solved
* *YL* liquid-phase mass fraction
* *YG* gas-phase mass fraction
* *YInt* value of mass fraction at the interface
*/

scalar Y[], YL[], YG[], YInt[];

/**
## User Data

Using this phase change model, the user should define the
following variables (SI units):

* *Dmix1* Diffusivity coefficient in the liquid phase
* *Dmix2* Diffusivity coefficient in the gas phase
* *YIntVal* Value of mass fraction at the interface on the
gas phase side. */

extern double Dmix1, Dmix2;
extern double YIntVal;

/**
Allocating diffusivity fields *Dmixf*, volume correction
*theta* and explicit source terms *sgS* and *slS* for
the diffusion equation. */

face vector Dmix1f[], Dmix2f[];
scalar thetacorr1[], thetacorr2[];
scalar sgS[], slS[];
scalar slimp[], sgimp[];

/**
## Defaults

In the defaults event we setup the tracer lists for the
advection of the mass fraction fields. */

event defaults (i = 0,last)
{
  YL.inverse = false;
  YG.inverse = true;

#ifdef CONSISTENTPHASE1
  fuext.tracers = list_append (fuext.tracers, YL);
#else
  fu.tracers = list_append (fu.tracers, YL);
#endif
#ifdef CONSISTENTPHASE2
  fuext.tracers = list_append (fuext.tracers, YG);
#else
  fu.tracers = list_append (fu.tracers, YG);
#endif
}

/**
## Init

In the init event, we avoid dumping all the fields that we
don't need to visualize. */

event init (i = 0)
{
  sgS.nodump = true;
  slS.nodump = true;
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
    sgS[] = 0.;
    slS[] = 0.;
    sgimp[] = 0.;
    slimp[] = 0.;
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
    YL[] = f[] > F_ERR ? YL[]/f[] : 0.;
    YG[] = ((1. - f[]) > F_ERR) ? YG[]/(1. - f[]) : 0.;
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
    mEvap[] = 0.; YInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      YInt[] = YIntVal;

      double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YIntVal, false);

#ifndef DIFFUSIVE
      mEvap[] = -rho2*Dmix2*gtrgrad/min(1. - YIntVal, 0.99);
#else
      mEvap[] = -rho2*Dmix2*gtrgrad;
#endif
    }
  }

  /**
  The source term for the diffusion equation of the
  species mass fractions in gas an liquid phase
  are computed here. */

  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fL, fsL), p;
      double alpha = plane_alpha (fL[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

#ifdef AXI
      sgS[] += -mEvap[]/rho2*area*(y + p.y*Delta)/(Delta*y)*cm[];
      slS[] +=  mEvap[]/rho1*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sgimp[] +=  mEvap[]/rho2*area*(y + p.y*Delta)/(Delta*y)*cm[];
      slimp[] += -mEvap[]/rho1*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      sgS[] += -mEvap[]/rho2*area/Delta*cm[];
      slS[] +=  mEvap[]/rho1*area/Delta*cm[];
      sgimp[] +=  mEvap[]/rho2*area/Delta*cm[];
      slimp[] += -mEvap[]/rho1*area/Delta*cm[];
#endif
    }
  }

  /**
  We restore the tracer form of the liquid and gas-phase
  mass fraction fields. */

  foreach() {
    YL[] *= f[]*(f[] > F_ERR);
    YG[] *= (1. - f[])*((1. - f[]) > F_ERR);
    Y[]  = YL[] + YG[];
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

We solve the diffusion equations for the mass fraction fields
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
    YL[] = fuext[] > F_ERR ? YL[]/fuext[] : 0.;
    YG[] = ((1. - fu[]) > F_ERR) ? YG[]/(1. - fu[]) : 0.;
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
    Dmix1f.x[] = Dmix1*fsL.x[]*fm.x[];
    Dmix2f.x[] = Dmix2*fsG.x[]*fm.x[];
  }

  foreach() {
    thetacorr1[] = cm[]*max(fL[], F_ERR);
    thetacorr2[] = cm[]*max(fG[], F_ERR);
    //slS[] = (f[] > F_ERR) ? slS[] : 0.;
    //sgS[] = (f[] > F_ERR) ? sgS[] : 0.;
    //slimp[] = (f[] > F_ERR) ? slimp[] : 0.;
    //sgimp[] = (f[] > F_ERR) ? sgimp[] : 0.;
  }

#ifndef SOLVE_LIQONLY
  diffusion (YG, dt, D=Dmix2f, r=sgS, beta=sgimp, theta=thetacorr2);
#endif
#ifndef SOLVE_GASONLY
  diffusion (YL, dt, D=Dmix1f, r=slS, beta=slimp, theta=thetacorr1);
#endif

  foreach() {

    /**
    Reconstruct mass fractions computing the mass of each
    chemical species in every cell. */

    double totmassliq = rho1*YL[]*fL[]*dv();
    YL[] = (totmassliq > 0.) ? rho1*YL[]*fL[]*dv() / totmassliq : 0.;

    double totmassgas = rho2*YG[]*fG[]*dv() + rho2*(1. - YG[])*fG[]*dv();
    YG[] = (totmassgas > 0.) ? rho2*YG[]*fG[]*dv() / totmassgas : 0.;

    YL[] *= fL[];
    YG[] *= (1. - fL[]);
  }

  /**
  We reconstruct the one-field mass fraction field summing
  the two fields $Y_L$ and $Y_G$ in tracer form. */

  foreach() {
    YL[] = clamp (YL[], 0., 1.);
    YG[] = clamp (YG[], 0., 1.);
    YL[] = (YL[] > F_ERR) ? YL[] : 0.;
    YG[] = (YG[] > F_ERR) ? YG[] : 0.;
    Y[] = YL[] + YG[];
  }
}

