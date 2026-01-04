/**
# Convergence of axisymmetric EHD stresses

We wish to test the accuracy of the EHD stresses in case of
axisymmetric coordinates. To do so, we create an "ad hoc" solution for
the electric poisson equation, 
$$ 
\frac {\partial_y (y \epsilon \partial_y \phi)}{y} + 
\partial_x(\epsilon \partial_x \phi) = - \rho_e
$$ 
with $\rho_e$ the electric charge density. If we look for a solution
of the form, 
$$ 
\phi (x, y) = y^a \cos x \sin y
$$ 
the charge density must then be (for $\epsilon = 1$), 
$$ 
\rho_e = -y^{-2 + a} \cos x \left( (1 + 2 a) y \cos y 
+ (a^2 - 2 y^2) \sin y \right) 
$$ 

The corresponding electrical stresses $\mathbf{F}$ are calculated from the
divergence of the Maxwell stress tensor, 
$$ 
\mathbf{F} = \nabla \cdot \mathbf{M} 
$$ 
Their components in the x- and y- directions are, 
$$
\mathbf{F}|_y = \frac{\partial_y (y M_{yy})}{y} + \partial_x M_{yx} -
\frac{M_{\theta \theta}}{y}
$$ 
and 
$$ \mathbf{F}|_x = \frac{\partial_y
(y \partial_y M_{xy})}{y} + \partial_x M_{xx}
$$ 
with 
$$
M_{x x} = \epsilon (E_x^2 - E_y^2)/2 
$$
$$ 
M_{y y} = \epsilon (E_y^2 - E_x^2)/2
$$ 
$$ 
M_{x y} = M_{y x} =\epsilon (E_y E_x) 
$$ 
and 
$$ 
M_{\theta \theta} = -\epsilon (E_y^2 + E_x^2)/2 
$$

$E_{x/y} = -\partial_{x/y} \phi$ are the axial and radial components 
of the electric field.

In this test we will set the charge density to check how the electric 
potential $\phi$ and the electrical stresses are recovered. */

#include "grid/multigrid.h"
#include "axi.h"
#include "run.h"
#include "poisson.h"

mgstats mg;
scalar phi[];
face vector epsilon[];
face vector a[], alpha[];
 
/**
We solve only for the case $a=1$. The case $a=0$ is ill-posed because
the forcing terms diverge near the axis of symmetry. */

#define PHI(x,y) y*cos(x)*sin(y)
#define RHOE(x,y) cos(x)*(-3*cos(y) + (2*y-1/y)*sin(y))

#define FY(x,y) (y == 0. ? 0. : sq(cos(x))/y*(sq(sin(y)) + \
					      y/2*(y + 5*y*cos(2*y)	\
						   - 2*(-2+sq(y))*sin(2*y))))

#define FX(x,y) cos(x)*sin(x)*sin(y)*(-3*y*cos(y)+(2*sq(y)-1)*sin(y))

phi[top] = dirichlet(0.);
phi[bottom] = dirichlet(0.);

int main()
{
  L0 = 2.*pi;
  periodic (right);
  for (N = 16; N <= 128; N *= 2)
    run();
}

event init (i = 0) {
  scalar rhoe[], rhs[];
  foreach() {
    rhoe[] = RHOE(x,y);
    rhs[] = -rhoe[]*cm[];
    phi[] = 0.;
  }

  foreach_face() {
    epsilon.x[] = fm.x[];
    alpha.x[] = fm.x[];
    a.x[] = 0.;
  }
  
  mg = poisson (phi, rhs, epsilon);
}

#include "ehd/stress.h"

event plot_err (i = 0)
{     
  face vector err[];

  foreach_face (x) {
    err.x[] = a.x[] - FX(x,y);
    printf ("x: %g %g %g %g %g  \n", x, y, a.x[], FX(x,y), err.x[]);
  }
     
  foreach_face (y) {
    err.y[] = a.y[] - FY(x,y);
    printf ("y: %g %g %g %g %g \n", x, y, a.y[], FY(x,y),err.y[]);
  }
  
  norm nx = normf (err.x), ny = normf (err.y);
     
  fprintf (stderr, "%d %g %g %g %g %g %g\n", N,
	   nx.avg, nx.rms, nx.max,
	   ny.avg, ny.rms, ny.max);
}

/**
As expected we get second-order convergence.

~~~gnuplot Convergence
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
fx(x)=a+b*x
fit fx(x) 'log' u (log($1)):(log($4)) via a,b
fx2(x)=a2+b2*x
fit fx2(x) 'log' u (log($1)):(log($2)) via a2,b2
fy(x)=c+d*x
fit fy(x) 'log' u (log($1)):(log($7)) via c,d
fy2(x)=c2+d2*x
fit fy2(x) 'log' u (log($1)):(log($5)) via c2,d2
set xlabel 'Resolution'
set ylabel 'Error'
set logscale
set xrange [8:256]
set cbrange [1:2]
set xtics 8,2,256
set grid ytics
set yrange [:]
plot 'log' u 1:4 t 'Fx:max', exp(fx(log(x))) t ftitle(a,b), \
     'log' u 1:2 t 'Fx:norm1', exp(fx2(log(x))) t ftitle(a2,b2), \
     'log' u 1:7 t 'Fy:max', exp(fy(log(x))) t ftitle(c,d),     \
     'log' u 1:5 t 'Fy:norm1', exp(fy2(log(x))) t ftitle(c2,d2)
~~~
*/
