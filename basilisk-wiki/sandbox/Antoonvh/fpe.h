/**  
# A Fokker-Planck equation solver

The Fokker-Plank equation (in 1D) reads, 

$$\frac{\partial \rho}{\partial t} = -\frac{\partial D_1 \rho}{\partial
x} + \frac{1}{2}\frac{\partial^2 D_2 \rho}{\partial x^2}$$

Which a candidate for a flux-divergence equation:

$$\frac{\partial \rho}{\partial t} = -\nabla \cdot J,$$

with,
$$J = D_1 \rho - \frac{1}{2} \frac{\partial D_2 \rho}{\partial x}$$

We asssume $D2$ to be a scalar field for now.
 */
#include "run.h"
scalar rho[], D2[];
face vector D1[];

void fpe_flux (face vector J, scalar rho,
	       (const) face vector D1, (const) scalar D2) {
  foreach_face() {
    J.x[] =  D1.x[]*(D1.x[] > 0 ? rho[-1] : rho[]) - //Fist order upwind
      (D2[]*rho[] - D2[-1]*rho[-1])/(2.*Delta);
    if (fabs(x - X0) < Delta/2)
      J.x[] = 0;
    else if (fabs(x - X0 - L0) < Delta/2.)
      J.x[] = 0;
  }
}

// Corresponding tendency

void fpe_tendency (scalar a, scalar da) {
  face vector J[];
  fpe_flux(J, a, D1, D2);
  foreach() {
    da[] = 0;
    foreach_dimension()
      da[] += -(J.x[1] - J.x[]);
    da[] /= Delta;
  }
}
/**
   We can use some integration scheme that uses these functions.
 */
#include "GLrk.h"  // Fix me: Not efficient
double fpe_TOL = 1e-8;
void wrap_fpe_tendency (scalar * al, scalar * dal) {
  scalar a, da;
  for (a, da in al, dal) //Fix me: Not efficient
    fpe_tendency(a, da);
}

event init (t = 0);

/**
The timestep is set by `DT` and the timing of other events.
 */

event timestep (i++) {
  dt = dtnext (DT);
}

event advance_fpe (i++, last) {
  A_Time_Step ({rho}, dt, wrap_fpe_tendency, fpe_TOL);
  //runge_kutta ({rho}, t, dt, wrap_fpe_tendency, 2);
}

