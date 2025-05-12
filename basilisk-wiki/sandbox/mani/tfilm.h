/**under develoment*/

/**
The approximate equations of motion valid in he limit of slender jets i.e. 
jets for which the radial dimensions are much smaller than the longi-tudinal dimensions will take the form
$$
\frac{\partial h}{\partial t} + v \frac{\partial h}{\partial z}  = 
-\frac{1}{2} h \frac{\partial v}{\partial z}  
$$
$$
\frac{\partial v}{\partial t} + v \frac{\partial v}{\partial z} = 
- \frac{\gamma}{\rho} \frac{\partial p}{\partial z} + 3\frac{\mu}{\rho*h^2} \frac{\partial}{\partial z}\left(h^2 \frac{\partial v}{\partial z} \right) 
$$ 
where $h(z,t)$ and $v(z,t)$ stand respectively for the radius of the
jet and the local velocity in a thin slice.  The ``pressure'' $p$
actually only amounts to the capillary contribution of pressure, and
is fully determined by the geometry of the jet:
$$
p(h) = \frac{1}{h \left(1 + h_z^2\right)^{1/2}} - \frac{h_{zz}}{\left(1+h_z^2\right)^{3/2}}.
$$
Note that here, the full expression for the curvature is used (see e.g
[Eggers and Dupont, J. Fluid Mech., 262
(1994)](http://dx.doi.org/10.1017/S0022112094000480) for a
discussion.)

## Conservative form

This set of equations can be rewritten in the following conservative form:
$$
\frac{\partial h^2}{\partial t} + \frac{\partial h^2 v}{\partial z}  = 0
$$
$$
\frac{\partial h^2 v}{\partial t} + \frac{\partial h^2 v^2}{\partial z} = 
\frac{\gamma}{\rho} \frac{\partial}{\partial z}\left(h^2 K\right) + 3\frac{\mu}{\rho} \frac{\partial}{\partial z}\left(h^2 \frac{\partial v}{\partial z} \right) 
$$
with $K$ denoting the following quantity:
$$
K = \frac{h_{zz}}{\left(1+h_z^2\right)^{3/2}} + \frac{1}{h \left(1 + h_z^2\right)^{1/2}}.
$$
(yes, this is now a `+' sign -- see [Li and Fontelos, Phys. Fluids,
15(4) (2003)](http://dx.doi.org/10.1063/1.1556291).)

This formulation reveals the mass per unit length $\rho h^2$ and the
momentum per unit length $\rho h^2 v$ of a thin slice as natural
variables of the problem. 

 
We will now solve the governing equations in the unstable regime, extract
the instability growth rate and compare it with the theoretical prediction.

We use a 1D Cartesian grid and the generic loop. */

#include "grid/cartesian1D.h"
#include "run.h"
/**
## Variables

We define the conservative variables as a *scalar* mass per unit length
and a *face vector* momentum per unit length. */

scalar h2[];
face vector h2u[];

/**Declare the properties density, surface tension coefficient, radius of the film, 
   dynamic viscosity, kinematic viscosity and define laplace number*/
 

double rho, mu, sigma, rad;
double gam = mu/rho;

/**#define LA sigma*rad/(rho*gamma*gamma)*/


/**
## Boundary conditions

By default symmetry conditions are applied on h2 and h2u 

*/


/**
## Explicit resolution

The integration is fully explicit. Note that the size of timestep is adjusted by 
computing the CFL and the distance to the next event. The $K$ term is computed by composing
stencil operations, and the convective flux is approximated in the Harlow-Welch style.

*/

event integration (i++) {
  double dt = DT;
  scalar u[];
  scalar h2uflux[];
  scalar d1u[];
  
  foreach_face() {
    u[] = 2.*h2u.x[]/(h2[] + h2[-1]);
    h2uflux[] = sq(h2u.x[1] + h2u.x[])/(4.*h2[]);
    double un = fabs(u[]);
    if (un > 0. && CFL*Delta*Delta/un < dt)
      dt = CFL*Delta*Delta/un;
  }
  dt = dtnext (dt);

  boundary({h2uflux});
  
  foreach_face() {
  d1u[]= (u[1] - u[])/Delta;
  }
  
  boundary({d1u});
  
  scalar h[];
  foreach()
    h[] = sqrt(h2[]);
  boundary ({h});

  scalar K[];
  foreach() {
    double hz = (h[1] - h[-1])/(2.*Delta);
    double hzz = (h[1] + h[-1] - 2.*h[])/sq(Delta);
    double ds = sqrt(1. + sq(hz));
    K[] = hzz/pow(ds,3) + 1./(h[]*ds);
  }
  boundary ({K});

  foreach_face()
    h2u.x[] += - dt*(h2uflux[] - h2uflux[-1])/Delta + dt*(sigma/rho)*(h2[]*K[] - h2[-1]*K[-1])/Delta + 1.5*dt*gam*((h2[1]+h2[])*d1u[] - (h2[1]+h2[-1])*d1u[-1])/Delta;

  scalar dh2 = K;
  foreach()
    dh2[] = - (h2u.x[1] - h2u.x[])/Delta;
  foreach()
    h2[] += dt*dh2[];

  boundary ({h2});
  boundary_flux ({h2u});
}

  
/** ##diagnostic equaiton for radial component of velocity */


