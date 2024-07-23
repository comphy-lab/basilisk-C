/**
# Two-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

#include "vof.h"

scalar f[], * interfaces = {f};

#include "two-phase-generic.h"

/**
## See also

* [Two-phase interfacial flows with coupled VOF and levelset](two-phase-clsvof.h)
* [Two-phase interfacial flows with levelset](two-phase-levelset.h)
*/
