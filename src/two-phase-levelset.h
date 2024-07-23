/**
# Two-phase interfacial flows with levelset

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a levelset
method. The signed distance to the interface is tracked by `d`. The
volume fraction field `f` is computed from the signed distance.

Note that the [coupled VOF and levelset solver](two-phase-clsvof.h)
should be prefered as it ensures much better mass conservation.

This solver can be combined with the [integral formulation of surface
tension](integral.h).

The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

#include "tracer.h"

scalar d[], f[], * tracers = {d};

#include "two-phase-generic.h"

void levelset_to_vof (scalar d, scalar f)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
  fractions (phi, f);
}

/**
The initial volume fraction is computed from the initial distance
field, which must be initialised by the user. */

event init (i = 0)
{
  levelset_to_vof (d, f);
}

/**
The distance function is reinitialised at each timestep. */

#include "redistance.h"

event properties (i++)
{
  redistance (d, imax = 3);
  levelset_to_vof (d, f);
}

/**
## See also

* [Two-phase interfacial flows with VOF](two-phase.h)
* [Two-phase interfacial flows with coupled levelset and VOF](two-phase-clsvof.h)
*/
