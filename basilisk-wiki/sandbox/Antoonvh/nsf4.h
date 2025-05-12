/**
# A 2D, face-averaged Navier-Stokes equations solver

For the solution to:

$$\frac{\partial \mathbf{u}}{\partial t} + \left(\mathbf{u} \cdot
\mathbf{\nabla}\right)\mathbf{u} = -\mathbf{\nabla} p + \nu \nabla^2
\mathbf{u} + \mathbf{a},$$

with the constraint that, 

$$\mathbf{\nabla} \cdot \mathbf{u} = 0.$$

## Method

The velocity components are represented discretely as face
*averages*. Their tendencies due to advection and diffusion are
computed at vertices wheareas the projection operators acts on the
face-averaged quantities.

It only works on periodic domains. :(
*/
#ifndef RKORDER
#define RKORDER (4)
#endif
#include "lsrkf.h"
#include "higher-order.h"
#include "my_vertex.h"

#ifdef PROJECT_2
#include "poisson.h"
#else
#include "../rajarshi/THESIS_CODES/Header_Files/poisson_O4.h"
#endif
/**
The global variables are,
*/
face vector u[];   // Face averaged values
scalar p[], p2[];  // Cell averaged values, not of the pressure!
(const) vector a;  // Vertex point values
(const) scalar nu; // Vertex point values
mgstats mgp, mgp2; // MG-solver statistics
/**
Some default settings that should work for most cases.
 */
event defaults (i = 0) {
#if TREE
  u.x.refine = refine_face_solenoidal; 
#ifndef PROJECT_2  
  TOLERANCE = 1e-3; 
  p.prolongation = refine_4th;
  p2.prolongation = refine_4th;
#endif
#endif
  foreach_dimension()
    u.x.prolongation = refine_face_4_x;
  if (a.x.i && !is_constant(a.x)) {
    foreach_dimension() {
      a.x.restriction = a.x.coarsen = restriction_vert;
      a.x.refine = a.x.prolongation = refine_vert5;
    }
  }
  CFL = 1.3;
  compact_iters = 5;
}

event init (t = 0);

event call_timestep (t = 0) {
  event ("timestep");
}
/**
## CFL timestepping

A stability criterion for the viscous term is neglected.
*/
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
An explicit 4th-order accurate second derivative scheme is used for
the viscous term.
 */
#define D2SDX2 (-(s[-2] + s[2])/12. + 4.*(s[1] + s[-1])/3. - 5.*s[]/2.)

/**
## The advection, viscous and acceleration terms
*/
void adv_diff (face vector du) {
  boundary ((scalar*){u});
  vector v[], dv[], dvx[], dvy[];
  scalar * vs = (scalar*){v, dv, dvx, dvy};
  for (scalar s in vs) {
    s.prolongation = refine_vert5;
    s.restriction = restriction_vert;
  }
  // Compute vertex values
  compact_face_av_to_vertex (u, v);
  // Compute derivatives
  compact_upwind ((scalar*){v}, {dvx, dvy}, v);
  // Apply a local projection and compute the tendency on vertices 
  foreach() {
    // Advection term
    double div = dvx.x[] + dvy.y[];
    dv.x[] = -(v.x[]*(dvx.x[] - div/2.) + v.y[]*dvx.y[]);
    dv.y[] = -(v.y[]*(dvy.y[] - div/2.) + v.x[]*dvy.x[]);
    // Viscous term
    if (nu.i) {
      foreach_dimension() {
	scalar s = v.x, ds = dv.x;
	foreach_dimension()
	  ds[] += nu[]*D2SDX2/sq(Delta);
      }
    }
    // Acceleration term
    if (a.x.i)
      foreach_dimension() 
	dv.x[] += a.x[];
  }
  // Compute the tendencies on faces
  boundary ((scalar*){dv});
  compact_vertex_to_face_av (dv, du);
}

/**
## Time integration

The well-known operator-splitting method is employed. Furthermore, we
keep track of the *worst* multigrid stratistics for all stages in
`mgp`.
*/

void Navier_Stokes (face vector u, face vector du) {
  adv_diff (du);
  boundary_flux({du});
  mgstats mgt = project (du, p, dt = dt);

  mgp.i      = max(mgp.i     , mgt.i);
  mgp.nrelax = max(mgp.nrelax, mgt.nrelax);
  mgp.resa   = max(mgp.resa  , mgt.resa);
  mgp.resb   = max(mgp.resb  , mgt.resb);
  mgp.sum    = max(mgp.sum   , mgt.sum);
}
/**
In order to prevent accumulation of the divergence' residuals, the
solution is projected after each itegration step.
 */
event advance (i++, last) {
  mgp = (mgstats){0}; // reset
  A_Time_Step (u, dt, Navier_Stokes); 
  boundary_flux ({u});
  mgp2 = project (u, p2, nrelax = 3); 
}

event adapt (i++, last);

/**
   ## Utilities
*/
#include "utils.h"

void vorticityf (face vector u, scalar omega) {
  foreach() {
    omega[] = 0;
    coord f = {1, -1};
    foreach_dimension()
      omega[] += f.x*(u.y[1,0] + u.y[1,1]
		      - u.y[-1,0] -u.y[-1,1]);
    omega[] /= (4*Delta);
  }
  boundary ({omega});
}

double Gauss6_x (double x, double y, double Delta, double (* myfun)(double x, double y)) {
  double w1 = 4./9., w2 = 5./18.;
  double yw = sqrt(3./5.)/2.*Delta;
  return w1*myfun (x, y) + w2*(myfun (x, y - yw) + myfun (x, y + yw));
}

double Gauss6_y (double x, double y, double Delta, double (* myfun)(double x, double y)) {
  double w1 = 4./9., w2 = 5./18.;
  double xw = sqrt(3./5.)/2.*Delta;
  return w1*myfun (x, y) + w2*(myfun (x - xw, y ) + myfun (x + xw, y));
}
/**
### A Wavelet-based grid-adaptation helper function

It can help reduce the likelyhood of many small/narrow
high-resolution islands. The power `p` can be used to optimize for a
certain error norm.
 */
#include "adapt_field.h"
astats adapt_flow (double ue, int maxlevel, double p) {
  boundary ((scalar*){u});
  scalar w[];
  vector wv[];
  foreach_dimension()
    wavelet (u.x, wv.x);
  for (int l = 1; l < depth(); l++) 
    foreach_coarse_level (l) { 
      double maxw = 0.;
      foreach_child()
	maxw = max(maxw, fabs(wv.x[]) + fabs(wv.y[]));
      w[] = maxw;
    }
  for (int l = depth(); l > 1; l--)  
    foreach_level(l)  // Gaussian kernel
      w[] = pow(Delta, p)*(2. *coarse(w,0,0,0)    +
			   1. *(coarse(w,1,0,0)   + coarse(w,0,1,0)   +
				coarse(w,-1,0,0)  + coarse(w,0,-1,0)) +
			   0.5*(coarse(w,1,1,0)   + coarse(w,1,-1,0)  +
				coarse(w,-1,-1,0) + coarse(w,-1,1,0)))/8.;
  return adapt_field (w, ue, ue/1.5, maxlevel);
}

/**
## Tests

* The 6th order Gaussian quadrature face diagnosis (Pass not shown). 
* The face-averages-to-vertices interpolation (Pass not shown)
* The vertices-to-face-averages interpolation (Pass not shown)
* The 4th order projection scheme (Pass not shown, see also Rajarshi Roy's thesis)
* The Vertex-based differencing schemes (Pass not shown)
* [Taylor-Green vortices (Pass)](tg4.c)
* [Adaptive stability (Pass)](ins4.c)

*/
