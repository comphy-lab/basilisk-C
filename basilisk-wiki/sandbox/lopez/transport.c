/**
# Concentration transport 

This test case intends to help in the understanding and evaluation of 
the transport of "VOF concentration" fields if the velocity field does 
not fulfill the incompressible velocity condition (the velocity 
is not divergence-free). This test case is complementary of [advect.c]().  */

#include "advection.h"
#include "src/vof.h"
#include "curvature.h"

/**
`ipos(t,a,xo)` compute the exact position of a planar interface
located initially at $x=xo$ if the velocity field is $\mathbf{u} = (x+a)
\mathbf{e_x}$. */

#define ipos(t, a, xo) (a*(exp(t)-1.) + xo*exp(t))

/**
The volume fraction is stored in scalar field `f` which is listed as
an *interface* for the VOF solver. We do not advect any tracer with
the default (diffusive) advection scheme of the advection solver. */

scalar f[], rho1[], rho2[];
scalar * interfaces = {f}, * tracers = NULL;
int MAXLEVEL = 5;

#define a 0.1
#define xo 0.

int main()
{

  /**
  We center the Box of length 2 on the origin.*/

  L0 = 2. [0];
  DT = HUGE [0];
  origin (-L0/2., -L0/2.);
  
  /**
  The scalar field `rho1` and `rho2` are "vof concentrations"
  associated with phase 1 determined by `f` obeying to the following
  equations,  
  $$ D_t \rho_1 = \partial_t \rho_1 + \mathbf{u} \cdot
  \nabla \rho_1 = 0 \quad \text{and} \quad \partial_t \rho_2 + \nabla
  \cdot (\mathbf{u} \rho_2) = 0 
  $$ 
  Note that both expression are not the same if 
  $\nabla \cdot \mathbf{u} \neq 0$. By default
  [vof.h](src/vof.h) solves always $D_t \rho = 0$ unless attribute
  *comp* were set to 1 (the default is zero).
  */

  f.tracers = {rho1, rho2};
  rho2.comp = 1;
  
  /**
  We then run the simulation for a certain refinement. */

  init_grid (1 << MAXLEVEL);
  run();
}

#define T 2.0

/**
Initially the fields`rho1` and `rho2` are homogeneous and equal to
1. The interface is located initially in $xo$, the velocity field is
unidirectional and linear in the x-direction and the timestep must be
set considering CFL. */

event init (i = 0)
{
  fraction (f, -x + xo);
  foreach() {
    rho1[] = f[];
    rho2[] = f[];
  }

  foreach_face(x)
    u.x[] = (x + a);
   dt = dtnext (timestep (u, DT));
}

/**
We compare the computed time evolution of the interface position with
the analytical solution. */

event interf_pos (t += 0.1) {
  scalar pos[];
  position (f, pos, {1,0});
  fprintf (stdout, "%g %g %g\n", t, statsf(pos).max, ipos(t, a, xo));
}

/**
At the final instant we compute the error made in the distribution in
both fields. */
  

event field (t = T) {

  scalar e[], ec[];
  foreach () {
    e[] = (f[] > 1e-12 ? rho1[]/f[]-1. : 0.);
    ec[] = (f[] > 1e-12 ? rho2[]/f[]-exp(-t) : 0.);
  }
  norm n = normf (e); 
  norm nc = normf (ec); 
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
  fprintf (stderr, "%d %g %g %g\n", N, nc.avg, nc.rms, nc.max);
  //dump ();
}
/**
~~~gnuplot Position of the interface.
set terminal @PNG enhanced size 640,640 font ",12"
set xlabel 't'
set ylabel 'interface position'
plot 'out' u 1:2  t 'vof.h', 'out' u 1:3 lt 2 w l t 'Analytical'
~~~
*/
