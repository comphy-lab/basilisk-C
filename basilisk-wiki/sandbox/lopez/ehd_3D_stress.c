/**
# Convergence of 3D EHD stresses

We wish to test the accuracy of the EHD stresses in 3D. 
To do so, we create an "ad hoc" solution for
the electric poisson equation, 
$$ 
{\partial_z (\epsilon \partial_z \phi)} + 
{\partial_y (\epsilon \partial_y \phi)} + 
\partial_x(\epsilon \partial_x \phi) = - \rho_e
$$ 
with $\rho_e$ the electric charge density. If we look for a solution
of the form, 
$$ 
\phi (x, y, z) = z \cos x \sin y \sin z
$$ 
the charge density must then be (for $\epsilon = 1$), 
$$ 
\rho_e = (3z \sin z - 2\cos z) \cos x \cos y 
$$ 

The corresponding electrical stresses $\mathbf{F}$ are calculated from the
divergence of the Maxwell stress tensor $\mathbf{M}$, 
$$ 
\mathbf{F} = \nabla \cdot \mathbf{M} 
$$ 

In this test we will set the charge density to check how the electric 
potential $\phi$ and the electrical stresses are recovered. */

#include "grid/multigrid3D.h"
#include "run.h"
#include "poisson.h"

mgstats mg;
scalar phi[];
const scalar rho[] = 1.;
scalar rhoe[], rhs[];
const face vector epsilon[] = {1.,1.,1.};
const face vector alpha[] = {1.,1.,1.};
face vector a[];
 
#define PHI(x,y,z) (z*cos(x)*sin(z)*sin(y))
#define RHOE(x,y,z) ((3*z*sin(z) - 2*cos(z))*cos(x)*sin(y))

#define FX(x,y,z) (z*cos(x)*sin(x)*sq(sin(y))*sin(z)*(-2*cos(z) + 3*z*sin(z)))
#define FY(x,y,z) ( - (3.*z*sin(z) - 2.*cos(z))*z*sq(cos(x))*sin(z)*sin(y)*cos(y))
#define FZ(x,y,z) ( - (3.*z*sin(z) - 2.*cos(z))*sq(cos(x))*sq(sin(y))*(sin(z) + z*cos(z)))

phi[front] = dirichlet(0.);
phi[back] = dirichlet(0.);

int main()
{
  L0 = 2.*pi;
  periodic (right);
  periodic (top);
  for (N = 16; N <= 128; N *= 2)
    run();
}

event init (i = 0) {
  foreach() {
    rhoe[] = RHOE(x,y,z);
    rhs[] = -rhoe[];
  }

  foreach_face()
    a.x[] = 0.;
  
  mg = poisson (phi, rhs, epsilon);
}

#include "src/stress.h"

event plot_err (i = 0)
{
  face vector err[];
  foreach_face (x)
    err.x[] = a.x[] - FX(x,y,z);
     
  foreach_face (y) 
    err.y[] = a.y[] - FY(x,y,z);

  foreach_face (z)
    err.z[] = a.z[] - FZ(x,y,z);

  norm nx = normf (err.x), ny = normf (err.y), nz = normf (err.z);
     
  fprintf (stderr, "%d %g %g %g %g %g %g %g %g %g\n", N,
  	   nx.avg, nx.rms, nx.max,
  	   ny.avg, ny.rms, ny.max,
  	   nz.avg, nz.rms, nz.max);
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
fz(x)=e+f*x
fit fz(x) 'log' u (log($1)):(log($10)) via e,f
fz2(x)=e2+f2*x
fit fz2(x) 'log' u (log($1)):(log($8)) via e2,f2
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
     'log' u 1:5 t 'Fy:norm1', exp(fy2(log(x))) t ftitle(c2,d2), \
     'log' u 1:10 t 'Fz:max', exp(fz(log(x))) t ftitle(e,f), \
     'log' u 1:8 t 'Fz:norm1', exp(fz2(log(x))) t ftitle(e2,f2)
~~~
*/
