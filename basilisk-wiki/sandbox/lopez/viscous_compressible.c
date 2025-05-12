/**
# Convergence of extra viscous terms appearing in compressible flows

We wish to test the accuracy of the discretisation of extra viscous
terms appearing in compressible fluids. To do so, we "manufacture" a
solution to the pure viscous diffusion problem
$$
  \partial_t (\mathbf{u}) = \nabla \cdot \tau + \mathbf{s}
$$
with
$$
  \tau = \mu (\nabla \mathbf{u} + \mathbf{u}^T) 
  + \lambda_v (\nabla \cdot \mathbf{u}) \mathbf{I}\, .
$$
and where the source term, $\mathbf{s} = s_x \mathbf{e}_x + s_y \mathbf{e}_y$, 
is to be determined. We look for solutions of the form
$$
u_x (x, y) = u_y (x, y) = y \cos x \sin y
$$
The viscous stress is then (for $\mu = 0$, $\lambda_v = 1$ and 
a axisymmetric domain), 
$$
\tau_{x x} = -s_x = - y \cos x \sin y  - (y \cos y  + 2 \sin y) \sin x \, ,
$$
$$
\tau_{y y} = -s_y =  \cos (y - x)  + 2  \cos (y + x)  - y \sin (y + x)
$$
and  $\tau_{x y} = 0$. For a 2D domain it would be,
$$
\tau_{x x} = -s_x = - y \cos x \sin y  - (y \cos y  + \sin y) \sin x \, ,
$$
and
$$
\tau_{y y} = -s_y =  \frac{1}{2}\left(\cos (y - x)  
+ 3  \cos (y + x)  - 2 y \sin (y + x) \right) \, .
$$

The game is then to solve the viscosity equation with the forcing
terms defined above and recover the (stationary) solution for the
velocity field. 

We use the axisymmetric metric, the viscous solver and the standard
time loop. */

#include "grid/multigrid.h"
#include "axi.h"
#define COMPRESSIBLE 1
#include "src/viscosity_compressible.h"
#include "run.h"

/**
We impose the correct boundary conditions for $u$ (the default
boundary conditions for $v$ are already correct). */

vector u[];

u.t[top] = dirichlet(0);
u.t[bottom] = dirichlet(0);

int main()
{
  L0 = 2.*pi;
  periodic (right);
  for (N = 16; N <= 128; N *= 2)
  run();
}

mgstats mg;

event init (i = 0) {
  foreach()
    u.x[] = u.y[]  = y*sin(y)*cos(x);
}

event integration (i++)
{
  double dt = 1.;

  /**
  This is the time integration loop. We first add the forcing term... */
  
  foreach() {
#if AXI
    u.x[] += dt*(y*cos(x)*sin(y) + y*cos(y)*sin(x) + 2.*sin(y)*sin(x));
    u.y[] += dt*(-cos(y-x) - 2.*cos(x+y) + y*sin(y+x));
#else
    u.x[] += dt*(y*cos(y)*sin(x) + (y*cos(x) + sin(x))*sin(y));
    u.y[] += dt*(sin(x)*(y*cos(y) + sin(y)) + cos(x)*(-2*cos(y) + y*sin(y)));
#endif
  }

  /**
  ...and then solve for viscosity. */
  
#if AXI
  mg = viscosity (u, mu = zerof, rho = cm, dt, lambdav = unityf);
#else
  mg = viscosity (u, mu = zerof, rho = unity, dt, lambdav = unityf);
#endif
}

event error (i = 20)
{

  /**
  Twenty timesteps are enough to converge toward the stationary
  solution. */
  
  scalar e[];
  foreach() {
    printf ("%g %g %g %g %g\n", x, y, u.x[], u.y[], y*sin(y)*cos(x));
    e[] = u.x[] - y*sin(y)*cos(x);
  }
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g %d %d\n",
	   N, n.avg, n.rms, n.max, mg.i, mg.nrelax);
}

/**
As expected we get second-order??? convergence. 

~~~gnuplot Convergence
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2
set xlabel 'Resolution'
set ylabel 'Error'
set logscale
set xrange [8:256]
set cbrange [1:2]
set xtics 8,2,256
set grid ytics
set yrange [1e-5:]
plot 'log' u 1:4 t 'max', exp(f(log(x))) t ftitle(a,b), \
     'log' u 1:2 t 'norm1', exp(f2(log(x))) t ftitle(a2,b2)
~~~

*/
