/**
# A 2D finite-diference co-located Navier-Stokes-equation solver

For the approximate solution to:

$$\frac{\partial \mathbf{u}}{\partial t} + \left(\mathbf{u} \cdot
\mathbf{\nabla}\right)\mathbf{u} = -\mathbf{\nabla} p + \nu \nabla^2
\mathbf{u},$$

with the constraint that, 

$$\mathbf{\nabla} \cdot \mathbf{u} = 0.$$

## Method

A Low-storage Runge-Kutta scheme is used to advance the centered
velocity field `u`. Two additional scalar fields are defined that
represent the pressure `p` and another field for extra projections
`p2`, which is not a pressure. Both are a solution to a Poisson
equation.
*/
#include "lsrk.h"
#include "poisson.h"
mgstats mgp, mgp2;
int it_project = 5; 
vector u[];
scalar p[], p2[];
double nu;

event defaults (i = 0) 
  CFL = 1.3;

mgstats projecter (vector v, scalar d, double Tol); // A prototype

event init (t = 0) 
  projecter (u, p2, TOLERANCE); //Project the user-initialized field

event call_timestep (t = 0) 
  event ("timestep");

/**
The timestep size is limited by the CFL condition *and* a Peclet
condition, the latter scales to be more stringent with decreasing
$\Delta$.

$$\frac{\nu \mathtt{dt}}{\Delta^2} \leq \mathcal{P}e$$
 */
double Pe = 0.4;
event timestep (i++, last) {
  double dtm = HUGE;
  foreach() {
    foreach_dimension() {
      if (fabs(u.x[]) > 0)
	if (Delta/fabs(u.x[]) < dtm)
	  dtm = Delta/fabs(u.x[]);
    }
  }
  dtm *= CFL;
  if (nu) {
    double Deltamin = L0/(1 << grid->maxdepth);
    dtm = min(dtm, Pe*sq(Deltamin)/nu);
  }
  dt = dtnext (min(DT, dtm));
}
/**
## The non-linear advection term

   For the advection term , 3rd-order accurate upwind differencing is
   used.
*/
#define ducx(uc) (u.x[] > 0 ? (2*uc[1]     + 3*uc[] -6*uc[-1]     + uc[-2])/6.     : \
		  (-2.*uc[-1]     - 3.*uc[] + 6.*uc[1]     - uc[2])/6.)
#define ducy(uc) (u.y[] > 0 ? (2*uc[0,1]   + 3*uc[] -6*uc[0,-1]   + uc[0,-2])/6.   : \
		  (-2.*uc[0,-1]   - 3.*uc[] + 6.*uc[0,1]   - uc[0,2])/6.)
#define ducz(uc) (u.z[] > 0 ? (2*uc[0,0,1] + 3*uc[] -6*uc[0,0,-1] + uc[0,0,-2])/6. : \
		  (-2.*uc[0,0,-1] - 3.*uc[] + 6.*uc[0,0,1] - uc[0,0,2])/6.)

void advect (vector du) {
  foreach()
    foreach_dimension()
      du.x[] = -(u.x[]*ducx(u.x) + u.y[]*ducy(u.x))/Delta;
}
/**
## Diffusion term

2nd-order accurate diffusion via the centered Lapacian.
 */
void diffuse (vector du) {
  foreach() 
    foreach_dimension() 
      du.x[] += nu*(u.x[0,1] + u.x[1,0] + u.x[0,-1] + u.x[-1,0] - 4.*u.x[])/sq(Delta);
}

/**
## Pressure-gradient term

The pressure is found as a solution to the Poisson equation:

$$\nabla^2 p = \left(\frac{\partial u_x}{\partial x}\right)^2 + 
2\frac{\partial u_x}{\partial y} \frac{\partial u_y}{\partial x} + 
\left(\frac{\partial u_y}{\partial y}\right)^2 $$

The right-hand side is computed using upwinding. The tendency due to
the pressure gradient is added to `du`.
*/
mgstats pressure_term (vector du, scalar p, double Tol) {
  scalar rhs[];
  foreach() {   
    rhs[] = -(sq(ducx(u.x)) +
	       2*ducy(u.x)*ducx(u.y) +
	      sq(ducy(u.y)))/sq(Delta);
  }
  mgstats mg = poisson (p, rhs, tolerance = Tol);
  foreach() 
    foreach_dimension()
      du.x[] -= (p[1] - p[-1])/(2*Delta);
  return mg;
}
/**
## The Runge-Kutta step

The functions above are called during the stages of the RK time
integrator.
 */
void NS_step (scalar * ul, scalar * dul) {
  vector du = {dul[0], dul[1]};
  advect (du);
  if (nu) 
    diffuse (du);
  mgp = pressure_term (du, p, TOLERANCE/sq(dt));
}
/**
## Time integration

Runge-Kutta Time integration is used. Because our methods are not
exactly solenoidal, an additional projection of the solution field
is performed every `it_project`-th iteration.
*/
mgstats projecter (vector v, scalar d, double Tol) {
  scalar div[];
  foreach() {
    div[] = 0;
    foreach_dimension()
      div[] += (v.x[1] - v.x[-1])/(2*Delta);
  }
  mgstats mg = poisson (d, div, tolerance = Tol);
  foreach() 
    foreach_dimension()
      v.x[] -= (d[1] - d[-1])/(2*Delta);
  boundary ((scalar*){v});
  return mg;
}

event advance (i++) {
  A_Time_Step ((scalar*){u}, dt, NS_step);
  boundary ((scalar*){u});
  if (i%it_project == 0)
    mgp2 = projecter (u, p2, TOLERANCE/sq(dt));
}

#if TREE
event adapt (i++, last);
#endif
/**
## Tests  

* [Dipole-Wall collision](tnsrk.c)
* [Merging of vortices](vortex.c)
* [Von Karman vortex street](karman.c)
* [Taylor-Green vortices](taylor-green.c)
*/
