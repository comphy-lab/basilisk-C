/**
# Two-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a Edge-Based Interface
Tracking (EBIT) method. 

The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. It's the same as that for VOF.*/

scalar f[], * interfaces = {f};
vector * u_ebit = {u};

#include "ebit-2d.h"

#include "two-phase-generic.h"

/**
## See also

* [Two-phase interfacial flows with VOF](two-phase.h)
* [Two-phase interfacial flows with levelset](two-phase-levelset.h)
* [Two-phase interfacial flows with coupled levelset and VOF](two-phase-clsvof.h)
*/
