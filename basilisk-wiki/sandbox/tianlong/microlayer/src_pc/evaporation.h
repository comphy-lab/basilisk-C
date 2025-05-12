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
extern double rho1, rho2;
extern face vector ufext;
extern scalar rho;

/**
We add a bool that allows the droplet volume changes to be disabled. */

bool is_shrinking;

/**
## Init event

We activate the shrinking and we set the dummy volume fractions to the
value of the vof field.
*/

event init (i = 0)
{
  is_shrinking = true;
  
  fu.nodump = true;
  fuext.nodump = true;

  foreach() {
    fu[] = f[];
    fuext[] = f[];
  }

  boundary((scalar *){fu,fuext,vpc});
}

/**
We add an event that can be overloaded by the user in order to
change some initialization parameters before starting the phase
change event (for example, initialize chemical species and temperature
to a specific initial solution).
*/

event end_init (i = 0);

/**
## Phase Change event

We computes the total vaporization rate per unit of surface
by summing all the contributions in *mEvapList*. We interpolate
the cell-centered normals and the mEvap field  in order to
compute a face change velocity *vpc*, defined at the edges.
*/

//use this to do some setup before phase change
event beforePhaseChange(i++);

event phasechange (i++)
{
  /**
  We restore the value of the *fu* and *fuext* fields
  at the beginning of the time step. */

  foreach() {
    fu[] = f[];
    fuext[] = f[];
  }
  boundary({fu,fuext});

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
    vpc.x[] = fm.x[]*(-mEvapf/rho2)*nf.x;
#else
    vpc.x[] = fm.x[]*(-mEvapf/rho1)*nf.x;
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
      coord m = normal (point, f);
      alpha = plane_alpha (f[], m);
      coord prel;
      segment = plane_area_center (m, alpha, &prel);
#ifdef AXI
      stefanflow[] = cm[]*mEvapTot[]*segment*(y + prel.y*Delta)*(1./rho2 - 1./rho1)/(Delta*y);
#else
      stefanflow[] = mEvapTot[]*segment*(1./rho2 - 1./rho1)/Delta*cm[];
#endif
    }
#endif
  }
}

/**
If *SHIFTING* is defined, we shift the expansion
source term. */

scalar stefanflowext[];

#ifdef SHIFTING
event shifting (i++) {
# ifdef EXT_STEFANFLOW
  foreach()
    stefanflowext[] = stefanflow[];
# endif
# ifdef SHIFT_TO_LIQ
  int dir = 1;
# else
  int dir = 0;
# endif
  shift_field (stefanflow, f, dir);
# ifdef EXT_STEFANFLOW
  shift_field (stefanflowext, f, 0);
# endif
}
#endif

/**
## VOF event

We overload the vof event in order to include the phase-change
velocity. Therefore, we perform the advection of the volume
fraction and we restore the original velocity field.
*/

static scalar * interfaces1 = NULL;

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

#if VELOCITY_JUMP
  extern face vector uf1;
  extern face vector uf2;
  extern vector n;
  extern scalar ls;
  extern scalar mEvapTot;
  extern scalar mEvapTotE;
  foreach_face()
  {
    double phif = 0.5 * (ls[] + ls[-1]);
    double mEvapf =  0.5 * (mEvapTotE[] + mEvapTotE[-1]);
    double nfx = 0.5 * (n.x[] + n.x[-1]);
    uf.x[] = phif > 0.0 ? uf1.x[] : uf2.x[];
    if(mEvapTot[] != 0.0 || mEvapTot[-1] != 0.0)
    {
      double rhof = phif > 0.0 ? rho1 : rho2;
      double jump = fm.x[] * mEvapf / rhof * nfx;
      uf.x[] += jump;
    }
  }
#endif
  //boundary((scalar *){uf});

  // It can be useful to compute the stability
  // conditions based on this modified velocity
  event ("stability");

  vof_advection ({f}, i);
  boundary ({f});

  /**
  We restore the value of the $\mathbf{uf}$ velocity field. */

  foreach_face()
    uf.x[] = uf_save.x[];
  boundary((scalar *){uf});

  /**
  We set the list of interfaces to NULL so that the default *vof()*
  event does nothing (otherwise we would transport $f$ twice). */

  interfaces1 = interfaces, interfaces = NULL;

#if VELOCITY_JUMP
  // here we compute the normal
  extern vector n;
  if (f.height.x.i > 0)
  {
    heights(f, f.height);
  }
  foreach ()
  {
    foreach_dimension()
        n.x[] = nodata;
    coord nf = normal(point, f);
    foreach_dimension()
        n.x[] = nf.x;
  }
#endif
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
  If the phase tracers are not advected with the same velocity
  of the volume fraction field, the advection is performed
  here exploiting the geometric advection procedure but
  using a dummy volume fraction field *fu* and *fuext*. */

  foreach_face() {
    uf_save.x[] = uf.x[];
    uf.x[] = ufext.x[];
  }
  boundary((scalar *){uf});

  /**
  We call the vof_advection function to transport
  the volume fraction *fuext* and associated tracers. */

  vof_advection ({fuext}, i);
  boundary({fuext});

  foreach_face()
    uf.x[] = uf_save.x[];
  boundary((scalar *){uf});

  /**
  We call the vof_advection function to transport
  the volume fraction *fu* and associated tracers. */

  vof_advection ({fu}, i);
  boundary({fu});

  /**
  We restore the original list of interfaces, which is
  required by [tension.h](src/tension.h). */

  interfaces = interfaces1;
}

