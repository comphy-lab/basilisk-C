/**
# Test of dissipation

We wish to test the accuracy of the discretisation of the viscous 
dissipation, $\Phi_v$,
$$
\Phi_v =  \tau : \nabla \mathbf{u}
$$
with
$$
  \tau = \mu (\nabla \mathbf{u} + \nabla \mathbf{u}^T) 
  + \lambda_v (\nabla \cdot \mathbf{u}) \mathbf{I}\, .
$$
With velocity profile of the form
$$
u_x (x, y) = u_y (x, y) = y \cos x \sin y
$$

The game is then compare the exact dissipation with the numerically 
computed one. */

#include "grid/multigrid.h"
#include "axi.h"
#define COMPRESSIBLE 1
#include "src/dissipation.h"
#include "run.h"

vector u[];
scalar dis[];

u.t[top] = dirichlet(0);
u.t[bottom] = dirichlet(0);

int main()
{
  L0 = 2.*pi;
  periodic (right);
  N = 32;
  for (N = 16; N <= 128; N *= 2)
  run();
}

event init (i = 0) {
  foreach()
    u.x[] = u.y[] = y*cos(x)*sin(y);

  dissipation (dis, u, fm, unityf);

  scalar e[];
  foreach() {
    double exact =  (2*y*sin(x)*sin(y)*(2*y*sin(x)*sin(y) - 
					cos(x)*(2.*y*cos(y) + 
						3.*sin(y))) +
		     sq(cos(x))*(4.*sq(y)*sq(cos(y)) +
				 9.*sq(sin(y)) + 5.*y*sin(2.*y)));
    fprintf (stdout,"%g %g %g %g\n", x, y, dis[], exact);
    e[] = dis[] - exact;
  }
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g \n",
	   N, n.avg, n.rms, n.max);
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
