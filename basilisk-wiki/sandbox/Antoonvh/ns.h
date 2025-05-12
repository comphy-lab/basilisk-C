/**
# A 2D finite-diference co-located Navier-Stokes-equation solver

A simple implemetation for the solution to:

$$\frac{\partial \mathbf{u}}{\partial t} + \left(\mathbf{u} \cdot
\mathbf{\nabla}\right)\mathbf{u} = -\mathbf{\nabla} p + \nu \nabla^2
\mathbf{u},$$

with the constraint that, 

$$\mathbf{\nabla} \cdot \mathbf{u} = 0.$$
 */
#include "run.h"

vector u[];
scalar p[];
double nu;

/**
The user is expected to initialize a flow. The timestep is limted
by the `CFL` condition.
 */
event init (t = 0);

event call_timestep (t = 0) {
  event ("timestep");
}

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
  dt = dtnext (min(DT, dtm));
}
/**
## The non-linear advection term

   For the advection term, 3rd order upwind differencing is used, with
   a slight upwind interpolation of the tendency.
 */
#define ducx(s, uc) (s > 0 ? (2*uc[1]     + 3*uc[] -6*uc[-1]     + uc[-2])/3.     : (-2.*uc[-1]     - 3.*uc[] + 6.*uc[1]     - uc[2])/3.)
#define ducy(s, uc) (s > 0 ? (2*uc[0,1]   + 3*uc[] -6*uc[0,-1]   + uc[0,-2])/3.   : (-2.*uc[0,-1]   - 3.*uc[] + 6.*uc[0,1]   - uc[0,2])/3.)
#define ducz(s, uc) (s > 0 ? (2*uc[0,0,1] + 3*uc[] -6*uc[0,0,-1] + uc[0,0,-2])/3. : (-2.*uc[0,0,-1] - 3.*uc[] + 6.*uc[0,0,1] - uc[0,0,2])/3.)
double UPW = 1.1; // Upwind interpolation parameter

void advect (double dt) {
  vector du[];
  foreach()
    foreach_dimension()
      du.x[] = (  u.x[]*(ducx(u.x[], u.x)) 
		+ u.y[]*(ducy(u.y[], u.x))
#if DIMENSION == 3
		+ u.z[]*(ducz(u.z[], u.x))
#endif
		)/(2.*Delta);
  boundary ((scalar*){du});

  foreach() {
    double CFLl = dt*(fabs(u.x[]) + fabs(u.y[]))/Delta;
    double iCFL = CFLl > 0 ? UPW/CFLl : HUGE;
    double xp = x - u.x[]*dt/(UPW + iCFL);
    double yp = y - u.y[]*dt/(UPW + iCFL);
    foreach_dimension()
      u.x[] -= dt*interpolate_linear(point, (struct _interpolate){du.x, xp, yp});
  }
  boundary ((scalar*){u});
}
/**
## Diffusion term

Both the forward and the backward Euler methods are available for the
diffusion term (2nd order in space).
 */
#if IMPLICIT_DIFFUSION
#include "diffusion.h"
#endif
void diffuse (double dt) {
#if !IMPLICIT_DIFFUSION
  vector du[];
  foreach() 
    foreach_dimension() 
      du.x[] = nu*(u.x[0,1] + u.x[1,0] + u.x[0,-1] + u.x[-1,0] - 4.*u.x[])/sq(Delta);
  foreach()
    foreach_dimension()
      u.x[] += dt*du.x[];
  boundary ((scalar*){u});
#else
  const face vector mu[] = {nu, nu};
  foreach_dimension()
    diffusion (u.x, dt, mu);
#endif
}
/**
## The pressure-Gradient term and enforcing the constraint

The projection step requires to solve a linear system. We use the
multigrid accelerated approach provided by `poisson.h`.
*/

#include "poisson.h"
mgstats mgp;
void projecter (double dt) {
  scalar rhs[];
  foreach() {   
    rhs[] = 0;
    foreach_dimension()
      rhs[] += (u.x[1] - u.x[-1])/(dt*2.*Delta);
  }
  mgp = poisson (p, rhs, tolerance = TOLERANCE/sq(dt));
  foreach() 
    foreach_dimension()
      u.x[] -= dt*(p[1] - p[-1])/(Delta*2.);
  boundary ((scalar*){u});
}
/**
## Operator splitting

This is not the Classical Chorin method as the viscous term is
computed from some intermediate field after advection.
 */

event advect_diffuse_and_project (i++) {
  advect (dt);
  if (nu) 
    diffuse (dt);
  projecter (dt);
}

#if TREE
event adapt (i++);
#endif
/**
## Test  

* [A non critical test](tns.c)
*/
