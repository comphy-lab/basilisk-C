/**
# A 2D, face-averaged Navier-Stokes equations solver

For the solution to:

$$\frac{\partial \mathbf{u}}{\partial t} + \left(\mathbf{u} \cdot
\mathbf{\nabla}\right)\mathbf{u} = -\mathbf{\nabla} p + \nu \nabla^2
\mathbf{u},$$

with the constraint that, 

$$\mathbf{\nabla} \cdot \mathbf{u} = 0.$$

## Method

The velocity components are represented discretely as face
*averages*. Their tendencies due to advection and diffusion are
computed at vertices wheareas the projection operator acts on the
prognostic face-averaged quantities.

## Options for the approximations

Options are:

* 3rd or 4th-order accurate time integration  
    - set `RKORDER`, default `RKORDER = 3` 
* 2nd-upwind, 3rd-upwind or 4th-order central differencing for the advection term 
    - swtich `ADV_4` or ADV_3, default is 2nd order
    - Switch `ADV_C` for combination: `max (ADV_4, ADV_3)`
* 2nd or 4th-order accurate projection (See Rajarshi Roy's work) 
    - switch `PROJECT_4`, the default is 2nd order
*/
#include "lsrkf.h"
#include "higher-order.h"
#include "my_vertex.h"

#ifdef PROJECT_4
#include "../rajarshi/THESIS_CODES/Header_Files/poisson_O4.h"
#else
#include "poisson.h"
#endif

face vector u[];
scalar p[], p2[];
mgstats mgp, mgp2;
double nu;

event defaults (i = 0) {
#if TREE
  u.x.refine = refine_face_solenoidal;
#ifdef PROJECT_4  //This improves convergence
  p.prolongation = refine_4th;
  p2.prolongation = refine_4th;
#endif
#endif
  foreach_dimension()
    u.x.prolongation = refine_face_4_x;
  CFL = 1.3;
  compact_iters = 5;
}

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

// 1st derivative schemes:
#define DSDX4 ((8*(s[1] - s[-1]) + s[-2] - s[2])/12)
#define DSDX3U (v.x[] > 0 ? \
		(2*s[1] + 3*s[] - 6*s[-1] + s[-2])/6. : \
		(-2*s[-1] - 3*s[] + 6*s[1] - s[2])/6)
#define DSDX4U (v.x[] > 0 ? \
	      (-s[-3] + 6*s[-2] - 18*s[-1] + 10*s[] + 3*s[1])/12.:	\
	      (-3*s[-1] - 10*s[] + 18*s[1] - 6*s[2] + s[3])/12.)
#define DSDX5U (v.x[] > 0 ?						\
		(-2.*s[-3] + 15*s[-2] - 60*s[-1] + 20*s[] + 30*s[1] - 3*s[2])/60.: \
		(3*s[-2] - 30*s[-1] -20*s[1] + 60*s[1] -15*s[2] + 2*s[3])/60.)
#define DSDX2 (v.x[] > 0 ? \
	       (3*s[] - 4*s[-1] + s[-2])/2. :	\
	       (-3*s[] + 4*s[1] - s[2])/2)
#ifdef ADV_4
#define DSDX (DSDX4)
#elif defined ADV_3
#define DSDX DSDX3U
#elif defined ADV_4U
#define DSDX (DSDX4U)
#elif defined ADV_C
#define DSDX (fabs(DSDX3U) > fabs(DSDX4) ? DSDX3U : DSDX4)
#else
#define DSDX DSDX2
#endif

// Diffusion
#define D2SDX2 (s[1] - 2*s[] + s[-1])


void adv_diff (face vector du) {
  boundary ((scalar*){u});
  vector v[], dv[];
  scalar * vs = (scalar*){v, dv};
  for (scalar s in vs) {
    s.prolongation = refine_vert5;
    s.restriction = restriction_vert;
  }
  compact_face_av_to_vertex (u, v);
  foreach() {
    foreach_dimension() {
      scalar s = v.x, ds = dv.x; // Don't rotate
      ds[] = 0;
      foreach_dimension() {
	ds[] -= v.x[]*DSDX/(Delta);
	if (nu)
	  ds[] += nu*D2SDX2/sq(Delta);
      }
    }
  }
  boundary ((scalar*){dv});
  compact_vertex_to_face_av (dv, du);
}
/**
## Time integration

Some operator-splitting method is employed. Furthermore, we keep track
of the *worst* multigrid stratistics for all stages in `mgp`.
*/

void Navier_Stokes (face vector u, face vector du) {
  adv_diff (du);
  mgstats mgt;
  mgp.i = 0;
  mgp.nrelax = 0;
  boundary_flux({du});
  mgt = project (du, p, dt = dt);
  mgp.i = max(mgp.i, mgt.i);
  mgp.nrelax = max(mgp.nrelax, mgt.nrelax);
}

event advance (i++, last) {
  A_Time_Step (u, dt, Navier_Stokes);
  boundary_flux ({u});
  mgp2 = project (u, p2, nrelax = 2); 
}

event adapt (i++, last);

/**
   ## Utilities
*/
#include "utils.h"

void vorticityf (face vector u, scalar omega) {
  foreach() {
    omega[] = 0;
    coord f = {1,-1};
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
## Tests

* The 6th order Gaussian quadrature face diagnosis (Pass not shown). 
* The face-averages-to-vertices interpolation (Pass not shown)
* The vertices-to-face-averages interpolation (Pass not shown)
* The 4th order projection scheme (Pass not shown, see also Rajarshi Roy's thesis)
* The Vertex-based differencing schemes (Pass not shown)
* [Taylor-Green vortices (Pass)](tg.c)
* [Adaptive stability (Pass)](ins.c)

*/
