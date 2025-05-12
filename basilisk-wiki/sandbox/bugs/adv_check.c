/**
#  Advection scheme mess up with tensors

*SP: This is not a bug. The default boundary conditions for scalar and
*tensor fields are different. They must be different because scalar,
*vector and tensor fields do not have the same symmetries. Setting the
*boundary conditions for the tensor components identically, see the
*code in main() below, gives the same results.*

Given the tracer fields $A$ and $B$ related by the function $B = \log
A$ ,for example.  Observe that if $A$ obeys to the advection equation
$B$ will do same,

$$
\quad \partial_t A + \mathbf{u} \cdot \nabla A = 0 \quad \Rightarrow \quad \partial_t B + \mathbf{u} \cdot \nabla B = 0 
$$

Therefore, the relationship $B = \log A$ should hold in time.

We will check that the advection scheme do not degenerate
the functional relationship as time goes by.  $A$ and/or $B$ 
could be scalars and/or components of a tensor. */

#define DT_MAX (0.001)
#define MU0 10. //  Viscosity 
#define uwall(x,t) (8*(1+tanh(8*(t-1/2)))*sq(x)*sq(1-x))

#include "navier-stokes/centered.h"


#include "tracer.h"
scalar A[], B[]; 
tensor T[];
scalar * tracers = {A, T.x.x, B};

int main()
{

  /**
  *SP: here we set the boundary conditions for T to that of a
  symmetric scalar field. We must use a temporary vector v because
  qcc is not clever enough to recognise T.x.t.* */
  
  vector v = T.x;
  v.t[top] = neumann(0);
  v.t[bottom] = neumann(0);
  v.t[left] = neumann(0);
  v.t[right] = neumann(0);
  
  DT = DT_MAX;
  N = 64;
  init_grid (N);
  const face vector mus[] = {MU0,MU0};
  mu = mus;
  run();
}

u.t[top] = dirichlet(uwall(x,t));
u.n[top] = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.n[bottom] = dirichlet(0);
u.t[left] = dirichlet(0);
u.n[left] = dirichlet(0);
u.t[right] = dirichlet(0);
u.n[right] = dirichlet(0);

event init (i = 0) 
{
  foreach() {
    u.x[] = 0.;
    A[] = (4+3*sin(2*x)*cos(2*y));
    B[] = log (A[]);
    T.x.x[] = log (A[]);
  }
  boundary ((scalar *){u, A, B, T.x.x});
}

static double energy()
{
  double se = 0.;
  if (u.x.face)
    foreach(reduction(+:se))
      se += (sq(u.x[] + u.x[1,0]) + sq(u.y[] + u.y[0,1]))/8.*sq(Delta);
  else // centered
    foreach(reduction(+:se))
      se += (sq(u.x[]) + sq(u.y[]))/2.*sq(Delta);
  return se;
}

event kinetic_energy (t += 0.01) 
{
  //  static FILE * fp = fopen ("kinetic", "w");
  fprintf (stderr, "%g %g\n", t, energy());
  //fflush(fp);
}

event profile (t = 1.8)
{

  FILE * fpp = fopen("yprof", "w");
  for (double y = 0; y <= 1.; y += 0.01)
    fprintf (fpp, "%g %g %g %g\n", y, 
	     log(interpolate (A, 0.5, y)), 
	     interpolate (B, 0.5, y), 
	     interpolate (T.x.x, 0.5, y)); 
  fclose (fpp);
  
  fpp = fopen("xprof", "w");
  for (double x = 0; x <= 1; x += 0.01)
    fprintf (fpp, "%g %g %g %g \n", x,
	     log(interpolate (A, x, 0.75)), 
	     interpolate (B, x, 0.75), 
  	     interpolate (T.x.x, x, 0.75));
  fclose (fpp);
}

/**

The scalars are correctly advected. The relationship $B = \log A$ holds at
instant $t=1.8$ . However, if the scalar is
a component of a tensor (the $xx$ component), the
relationship does not stand up, i.e $\mathbf{T}_{xx} \neq \log A$. It is
even worst if the component were $xy$-component.

~~~gnuplot x profile at time t=1.8 
plot 'xprof' u 1:2 w l t 'reference scalar', 'xprof' u 1:3 t 'scalar', 'xprof' u 1:4 t 'tensor'
~~~
*/
