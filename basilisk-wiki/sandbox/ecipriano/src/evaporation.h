/**
# Main Evaporation Module

Include this file every time you need to perform simulations
with phase change. It contains the declaration of fields
that are useful to every phase change model, it computes the
total evaporation rate and the phase change velocity, it performs
the stefan flow shifting, and it manages the advection of volume
fraction and tracers.
This file must be used in combination with a phase change model (i.e.
[multicomponent.h](multicomponent.h), 
[temperature-gradient.h](temperature-gradient.h), 
[species-gradient.h](species-gradient.h), 
[fixedflux.h](fixedflux.h)). */

#include "common-evaporation.h"
#include "fracface.h"
#include "diffusion.h"

/**
## Setup tracers advection

* *CONSISTENTPHASE1* (default): the liquid phase tracers are
advected using the velocity $\mathbf{u}^E$, associated
with the list fuext.tracers, while the gas phase tracers
are advected using $\mathbf{u}$ associated with the list fu.tracers.
* *CONSISTENTPHASE2*: the gas phase tracers are
advected using the velocity $\mathbf{u}^E$, associated
with the list fuext.tracers, while the liquid phase tracers
are advected using $\mathbf{u}$ using the list fu.tracers.
* *DIFFUSIVE*: the Stefan flow is neglected and no interface
velocity jump is present. Therefore, both the liquid and the
gas phase tracers are transported using $\mathbf{u}^E$ which
will be equal to $\mathbf{u}$, using the list fuext.tracers.
*/

#ifndef DIFFUSIVE
# ifndef CONSISTENTPHASE2
#   ifndef NOCONSISTENTPHASE
#     define CONSISTENTPHASE1
#   endif
# endif
#else
# define CONSISTENTPHASE1
# define CONSISTENTPHASE2
#endif

/**
## Setup Stefan flow

* *SHIFTING* (default): the volume expansion term is shifted.
* *NOSHIFTING*: do not apply shifting procedure.
* *SHIFT_TO_LIQ* (default): shifting toward the liquid phase (i.e. droplet evaporation).
* *SHIFT_TO_GAS*: shifting toward the gas phase (i.e. boiling problems).
*/

#define SHIFTING
#ifndef SHIFT_TO_GAS
# define SHIFT_TO_LIQ
#endif

#ifdef NOSHIFTING
# undef SHIFTING
#endif

/**
## Fields

We initialize fields useful to any phase change model:

* *mEvapTot*: total vaporization rate per unit of interface surface [kg/m2/s]
* *stefanflow*: expansion term used in the projection step of
[navier-stokes/centered-evaporation.h](navier-stokes/centered-evaporation.h)
* *vpc*: interface regression velocity [m/s]
* *fu*: dummy volume fraction used for the advection of scalar fields
using the directionally-split procedure, with the velocity $\mathbf{u}$
* *fuext*: dummy volume fraction used for the advection of scalar fields
using the directionally-split procedure, with the velocity $\mathbf{u}^E$
* *uf_save*: helper field used to store and recover a velocity before and
afer calling *vof_advection*
*/

scalar mEvapTot[];
scalar stefanflow[];
face vector vpc[];
scalar fu[], fuext[];
face vector uf_save[];

/**
## Extern fields

We initialize fields that must be defined externally. In particular,
*mEvapList* should be declared in the specific phase change model, and it
must contains a list of every contribution to the phase change (for example,
the vaporization rate of each chemical species).
*/

extern scalar * mEvapList;
extern face vector ufext;

/**
We add a bool that allows the droplet volume changes to be disabled. */

bool is_shrinking;
double smoothing_expansion = 0.;

/**
## Set Variable-Properties fields for compatibility with constant-properties
simulations. */

#ifndef VARPROP
//(const) scalar rho1v = unity, rho2v = unity;
//(const) scalar mu1v = zeroc, mu2v = zeroc;

// The (const) scalar approach gives problems when calling those variables
// as external by modules included before evaporation.h
scalar rho1v[], rho2v[];
scalar mu1v[], mu2v[];
#endif

/**
## Init event

We activate the shrinking and we set the dummy volume fractions to the
value of the vof field.
*/

event init (i = 0)
{
  is_shrinking = true;

  foreach() {
    fu[] = f[];
    fuext[] = f[];
  }

  scalar * interfaces2 = {fu,fuext};
#if TREE
  for (scalar c in interfaces2) {
    c.refine = c.prolongation = fraction_refine;
    c.dirty = true;
    scalar * tracers = c.tracers;
    for (scalar t in tracers) {
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = vof_concentration_refine;
      t.dirty = true;
      t.c = c;
    }
  }
#endif
  for (scalar c in interfaces2) {
    scalar * tracers = c.tracers;
    for (scalar t in tracers)
      t.depends = list_add (t.depends, c);
  }

#ifndef VARPROP
  //const scalar rho1vv[] = rho1, rho2vv[] = rho2;
  //rho1v = rho1vv, rho2v = rho2vv;

  //const scalar mu1vv[] = mu1, mu2vv[] = mu2;
  //mu1v = mu1vv, mu2v = mu2vv;

  foreach() {
    rho1v[] = rho1;
    rho2v[] = rho2;
    mu1v[] = mu1;
    mu2v[] = mu2;
  }
#endif
}

/**
We add an event that can be overloaded by the user in order to
change some initialization parameters before starting the phase
change event (for example, initialize chemical species and temperature
to a specific initial solution).
*/

event end_init (i = 0);

/**
## Reset Source Terms

We impose the position of the even which resets the source terms
for the governing equations solved in the specific phase change
models. */

event reset_sources (i++);

/**
## Chemical Reaction Step

If combustion chemistry is present we solve it here. */

event chemistry (i++);

/**
## Phase Change event

We computes the total vaporization rate per unit of surface
by summing all the contributions in *mEvapList*. We interpolate
the cell-centered normals and the mEvap field  in order to
compute a face change velocity *vpc*, defined at the edges.
*/

event phasechange (i++)
{
  /**
  We restore the value of the *fu* and *fuext* fields
  at the beginning of the time step. */

  foreach() {
    fu[] = f[];
    fuext[] = f[];
  }

  /**
  We compute the total vaporization mass-flowrate *mEvapTot*. */

  foreach() {
    mEvapTot[] = 0.;
    for (scalar mEvap in mEvapList)
      mEvapTot[] += mEvap[];
  }

  /**
  We update phase-change regression velocity *vpc*. */

  foreach_face() {
    vpc.x[] = 0.;
    coord nf;
    foreach_dimension()
      nf.x = 0.;
    double mEvapf = 0.;
    if (f[] > 0. && f[] < 1.) {
      nf = normal (point, f);
      mEvapf = mEvapTot[];
      if (f[-1] > 0. && f[-1] < 1.) {
        coord np = normal (neighborp(-1), f);
        foreach_dimension()
          nf.x = (nf.x + np.x)/2.;
        mEvapf = (mEvapf + mEvapTot[-1])/2.;
      }
    }
    else if (f[-1] > 0. && f[-1] < 1.) {
      mEvapf = mEvapTot[-1];
      nf = normal (neighborp(-1), f);
    }
#ifdef BYRHOGAS
    vpc.x[] = (rho2v[] > 0.) ? fm.x[]*(-mEvapf/rho2v[])*nf.x : 0.;
#else
    vpc.x[] = (rho1v[] > 0.) ? fm.x[]*(-mEvapf/rho1v[])*nf.x : 0.;
#endif
  }

  /**
  We compute the expansion source term for the
  continuity equation *stefanflow*. */

  foreach() {
    stefanflow[] = 0.;
#ifndef DIFFUSIVE
    double alpha = 0.;
    double segment = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord m = mycs (point, f);
      alpha = plane_alpha (f[], m);
      coord prel;
      segment = plane_area_center (m, alpha, &prel);
#ifdef AXI
      stefanflow[] = (rho1v[] > 0. && rho2v[] > 0.) ?
        cm[]*mEvapTot[]*segment*(y + prel.y*Delta)*(1./rho2v[] - 1./rho1v[])/(Delta*y)
        : 0.;
#else
      stefanflow[] = (rho1v[] > 0. && rho2v[] > 0.) ?
        mEvapTot[]*segment*(1./rho2v[] - 1./rho1v[])/Delta*cm[]
        : 0.;
#endif
    }
#endif
  }
}

/**
If *SHIFTING* is defined, we shift the expansion
source term. */

event shifting (i++) {

  /**
  We perform shifting operations of the expansion
  term, in order to remove it from the interfacial
  cells, but still keeping its effect on the velocity
  field. */

#ifdef SHIFTING
# ifdef SHIFT_TO_LIQ
  shift_field (stefanflow, f, 1);
# else
  shift_field (stefanflow, f, 0);
# endif
#endif

  /**
  To avoid an expansion term localized on a single layer
  of cells, we solve a diffusion equation which smoothen
  the expansion, and we solve the Projection step using
  this modified term. The volume integral of the expansion
  during this step must remain constant. */

  if (smoothing_expansion) {
    face vector eps[];
    foreach_face()
      eps.x[] = smoothing_expansion*fm.x[];

    scalar thetacorr[];
    foreach()
      thetacorr[] = cm[];

    diffusion (stefanflow, dt, D=eps);
  }
}

/**
## VOF event

We overload the vof event in order to include the phase-change
velocity. Therefore, we perform the advection of the volume
fraction and we restore the original velocity field.
*/

event vof (i++)
{
  foreach_face() {
    uf_save.x[] = uf.x[];
#ifndef INT_USE_UF
    uf.x[] = ufext.x[];
#endif
#ifndef CONSTANT_VOLUME
    if (is_shrinking)
      uf.x[] -= vpc.x[];
#endif
  }
}

/**
## Tracer Advection event

We use the directionally-split advection procedure in order
to resolve the convective transport of the scalar fields,
through the use of the dummy tracers *fu* and *fuext* which
exploit the velocity $\mathbf{u}_f$ and $\mathbf{u}_f^E$
respectively.
*/

event tracer_advection (i++)
{
  /**
  We restore the value of the $\mathbf{uf}$ velocity field. */

  foreach_face()
    uf.x[] = uf_save.x[];

  /**
  If the phase tracers are not advected with the same velocity
  of the volume fraction field, the advection is performed
  here exploiting the geometric advection procedure but
  using a dummy volume fraction field *fu* and *fuext*. */

  foreach_face() {
    uf_save.x[] = uf.x[];
    uf.x[] = ufext.x[];
  }

  /**
  We call the vof_advection function to transport
  the volume fraction *fuext* and associated tracers. */

  vof_advection ({fuext}, i);

#ifdef VELOCITY_JUMP
  foreach_face()
# ifdef BOILING_SETUP
    uf.x[] = uf1.x[];
# else
    uf.x[] = uf2.x[];
# endif
#else
  foreach_face()
    uf.x[] = uf_save.x[];
#endif

  /**
  We call the vof_advection function to transport
  the volume fraction *fu* and associated tracers. */

  vof_advection ({fu}, i);

#ifdef VELOCITY_JUMP
  foreach_face()
    uf.x[] = uf_save.x[];
#endif

}

