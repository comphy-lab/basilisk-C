/**
# Two-phase interfacial flows with coupled VOF and levelset

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The signed distance field `d` is advected as a tracer and is
relaxed toward the VOF-defined interface.

This coupling ensures a mass conservation at least as good as that of
the [pure VOF solver](two-phase.h).

This solver can be combined with the [integral formulation of surface
tension](integral.h).

The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

#include "vof.h"
#include "tracer.h"

scalar d[], f[], * interfaces = {f}, * tracers = {d};

#include "two-phase-generic.h"

/**
The initial volume fraction is computed from the initial distance
field, which must be initialised by the user. */

event init (i = 0)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
  fractions (phi, f);
}

/**
The distance function is reinitialised at each timestep. */

#include "redistance.h"

event properties (i++)
{

  /**
  In interfacial cells, the signed distance is obtained directly from
  the VOF reconstruction of the interface. This distance is combined
  with the existing distance using a small weight, thus ensuring
  exponential time relaxation of the signed distance toward its VOF
  value. */

  double weight = 0.1;
  foreach()
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f);
      normalize (&n);
      double alpha = plane_alpha (f[], n);
      d[] = (1. - weight)*d[] + weight*Delta*alpha;
    }

  /**
  The redistancing operation itself is quite expensive. */
  
  redistance (d, imax = 3);
}

/**
## See also

* [Two-phase interfacial flows with VOF](two-phase.h)
* [Two-phase interfacial flows with levelset](two-phase-levelset.h)
*/
