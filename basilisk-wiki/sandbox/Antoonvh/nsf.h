/**
# A 2D, finite-difference, staggered, Navier-Stokes-equation solver

A simple implemetation for the solution to:

$$\frac{\partial \mathbf{u}}{\partial t} + \left(\mathbf{u} \cdot
\mathbf{\nabla}\right)\mathbf{u} = -\mathbf{\nabla} p + \nu \nabla^2
\mathbf{u},$$

with the constraint that, 

$$\mathbf{\nabla} \cdot \mathbf{u} = 0.$$
 */
#include "run.h"
face vector u[];
scalar p[];
double nu;

u.n[left] = 0;
u.n[right] = 0;
u.n[bottom] = 0;
u.n[top] = 0;

/**
The user is expected to initialize a flow. The timestep is limited
by the `CFL` condition.
 */
event defaults (i = 0)
  CFL = .9;
		      
event init (t = 0);

event call_timestep (t = 0) {
  event ("timestep");
}

event timestep (i++, last) {
  double dtm = HUGE;
  foreach_face() {
    if (fabs(u.x[]) > 0)
      if (Delta/fabs(u.x[]) < dtm)
	dtm = Delta/fabs(u.x[]);
  }
  dtm *= CFL;
  dt = dtnext (min(DT, dtm));
}
/**
A function that applies the advection and diffusion operator for a
timestep `dt`:
*/
void adv_diff (double dt) {
  /**
   For the advection term, a centered tendency field is
   computed. Special care is taken to consider the staggering of the
   velocity-component fields.
  */
  vector du[];
  boundary ((scalar*){u}); 
  foreach() {
    double uc   = (u.x[] + u.x[1])/2;
    double vc   = (u.y[] + u.y[0,1])/2;
    double dudx = (u.x[1] - u.x[])/(Delta);
    double dudy = ((u.x[0,1] + u.x[1,1]) - (u.x[0,-1] + u.x[1,-1]))/(4*Delta);
    double dvdx = ((u.y[1,0] + u.y[1,1]) - (u.y[-1,0] + u.y[-1,1]))/(4*Delta);
    double dvdy = (u.y[0,1] - u.y[])/(Delta);
    du.x[] = uc*dudx + vc*dudy;
    du.y[] = uc*dvdx + vc*dvdy;
  }
  boundary ((scalar*){du});
  /**
`u` is first updated with the possible diffusion term. It is added in
the forward direction using the 5-point Laplacian.
   */
  if (nu) 
    foreach_face() 
      u.x[] += dt*nu*(u.x[0,1] + u.x[1,0] + u.x[0,-1] + u.x[-1,0]
		      - 4.*u.x[])/sq(Delta);
  /**
     The advection term (i.e. stored in `du`) is upwinded using interpolation. 
   */
  foreach_face(x) {
    double xp = x - u.x[]*dt/2.;
    double yp = y - (u.y[] + u.y[-1] + u.y[0,1] + u.y[-1,1])*dt/8.;
    u.x[] -= dt*interpolate_linear (point,
				    (struct _interpolate){du.x, xp, yp});
  }
  foreach_face(y) {
    double yp = y - u.y[]*dt/2.;
    double xp = x - (u.x[] + u.x[1] + u.x[0,-1] + u.x[1,-1])*dt/8.;
    u.y[] -= dt*interpolate_linear (point,
				    (struct _interpolate){du.y, xp, yp});
  }
}
/**
   The projection step requires the solution to a linear system. We
   use the multigrid-accelerated approach provided via `poisson.h`.
   The residual is weighted with $\mathtt{dt}^{-1}$, as suggested
   there.
 */
#include "poisson.h"
void projecter (double dt) {
  boundary ((scalar*){u});
  scalar div[];
  foreach() {   
    div[] = 0;
    foreach_dimension()
      div[] += u.x[1] - u.x[];
    div[] /= dt*Delta;
  }
  poisson (p, div);
  foreach_face() 
    u.x[] -= dt*(p[] - p[-1])/Delta;
}
/**
## Time integration

Chorin's classical operator-splitting method is employed. 
 */
event Navier_Stokes (i++) {
  adv_diff (dt);
  projecter (dt);
}

event adapt (i++, last);

/**
## Test  

* [A non critical test: Dipole-wall collision](tnsfh.c)
* [A Dipole-wall collision on a quadtree](tnsfa.c)
*/
