/**
# Convergence of axisymmetric viscous terms

We wish to test the accuracy of the discretisation of viscous terms in
axisymmetric coordinates. To do so, we "manufacture" a solution to the
pure viscous diffusion problem
$$
  \partial_t (yu) = \partial_x (y \tau_{x x}) + \partial_y (y \tau_{x y})
  + s_x
$$
$$
  \partial_t (yv) = \partial_x (y \tau_{y x}) + \partial_y (y \tau_{y y})
  + s_y - 2 \frac{v}{y}
$$
with
$$
  \tau_{x x} = 2\partial_x u
$$
$$
  \tau_{x y} = \tau_{y x} = (\partial_y u + \partial_x v)
$$
$$
  \tau_{y y} = 2\partial_y v
$$
and where $s_x$ and $s_y$ are source terms to be determined. We look
for solutions of the form
$$
u (x, y) = v (x, y) = y^a \cos x \sin y
$$
The viscous stress tensor is then (for $\mu = 1$)
$$
  \tau_{x x} = - 2 y^a \sin x \sin y
$$
$$
  \tau_{x y} = \tau_{y x} = \cos x (y^a \cos y + ay^{a - 1} \sin y) - y^a
  \sin x \sin y
$$
$$
  \tau_{y y} = 2 \cos x (y^a \cos y + ay^{a - 1} \sin y)
$$
and the viscosity equation gives
$$
  s_x = y^a  [3 y \cos x \sin y - (2 a + 1) \cos x \cos y - a^2 y^{- 1}
  \cos x \sin y + y \sin x \cos y + (a + 1) \sin x \sin y]
$$
$$
  s_y = y^a [y \sin x \cos y + a \sin x \sin y + 3 y \cos x \sin y - 2 (2
  a + 1) \cos x \cos y + 2 (1 - a^2) y^{- 1} \cos x \sin y]
$$
The game is then to solve the viscosity equation with the forcing
terms defined above and recover the (stationary) solution for the
velocity field. 

We use the axisymmetric metric, the viscous solver and the standard
time loop. */

#include "grid/multigrid.h"
#include "axi.h"
#include "viscosity.h"
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

  /**
  When using single precision (32-bits floats) convergence cannot be
  achieved to the default 1e-3. */
  
#if SINGLE_PRECISION
  TOLERANCE = 1e-2;
#endif

  for (N = 16; N <= 128; N *= 2)
    run();
}

/**
We solve only for the case $a=1$. The case $a=0$ is ill-posed because
the forcing terms diverge near the axis of symmetry. */

double a = 1.;
mgstats mg;

event init (i = 0) {
  foreach()
    u.x[] = u.y[] = 0.;
}

event integration (i++)
{
  double dt = 1.;

  /**
  This is the time integration loop. We first add the forcing
  term... */
  
  foreach() {
    u.x[] += dt*pow(y, a - 1.)*(3.*y*cos(x)*sin(y) - (2.*a + 1.)*cos(x)*cos(y)
				- sq(a)*cos(x)*sin(y)/y	
				+ y*sin(x)*cos(y) + (a + 1.)*sin(x)*sin(y));
    u.y[] += dt*pow(y, a - 1.)*(y*sin(x)*cos(y) + a*sin(x)*sin(y)
				+ 3.*y*cos(x)*sin(y)
				- 2.*(2.*a + 1.)*cos(x)*cos(y) +
				2.*(1. - sq(a))*cos(x)*sin(y)/y);
  }

  /**
  ...and then solve for viscosity. */
  
  mg = viscosity (u, fm, cm, dt);
}

event error (i = 20)
{

  /**
  Twenty timesteps are enough to converge toward the stationary
  solution. */
  
  scalar e[];
  foreach (serial) {
    printf ("%g %g %g %g %g\n", x, y, u.x[], u.y[], pow(y,a)*cos(x)*sin(y));
    e[] = u.x[] - pow(y,a)*cos(x)*sin(y);
  }
  norm n = normf (e);
  fprintf (stderr, "%d %.3g %.3g %.3g %d %d\n",
	   N, n.avg, n.rms, n.max, mg.i, mg.nrelax);
}

/**
As expected we get second-order convergence.

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
