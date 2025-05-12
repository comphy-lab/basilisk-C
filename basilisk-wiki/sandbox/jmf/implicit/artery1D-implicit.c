/**
# 1D arterial flow - implicit */

#include "grid/multigrid1D.h"
#include "artery-implicit.h"

/**
A 1D model for arterial flows can be derived from the Navier-Stokes
equations, in terms of the cross sectional area $A$ and flow rate $Q$,
we have

$$	 \partial_t A +\partial_x Q  = 0$$
$$  	 \partial_t Q +\partial_x (Q^2/A) = - A \partial_x p/\rho - f_r  $$
  	 
where $p(A)$ models the wall properties of the arteries, $\rho$ is the
blood density and $f_r$ stands by the wall shear stress. For a simple
linear wall relation, $p =2 K (A - 1)$ with $K$ constant we have the
same example than the explicit code.

Jose, y'a pas un bug? 

c'est $p = K (\sqrt{A} - 1)$ et ensuite la pression est une fausse
pression $\Pi$
$$  	 \partial_t Q +\partial_x (Q^2/A) = -  \partial_x \Pi - f_r  $$
donc $\Pi= K ( A^{3/2} - 1)/3$ car $\partial_x \Pi = \frac{K}{2}
A^{1/2} A_x = A \partial_x ( K A^{1/2})= A \partial_x ( p)$
*/

/**
The other parameters are specific to the example.
*/

double e1 ;
double e2 ;
double omega;
double Amp;


/**
## Functions

we define the functions  */

void update_source (scalar * current, scalar * updates)
{
  scalar a = current[0], q = current[1], da= updates[0], dq = updates[1];

  /**
  We add the source term for a and q.*/

  foreach() {
    da[] = 0.0 ;
    dq[] = -e2 * q[]/a[];
    // dq[] = -e2 * (q[] + q[-1])/(a[] + a[-1]);
  }
  boundary({da,dq});
 
}

/**
Functions pressure and derivative */

double pressurevalue(double h)
{
  return 2. * e1 * (pow(h,1.) - 1.);
}

double derivativevalue(double h)
{
  return 2. * e1;
}

/**
## Boundary conditions

We impose a sinusoidal flux Q(t) at the left of the domain. */

q[left] = dirichlet(Amp*sin(2.*pi*omega*t));
q[right] = neumann(0.0);

/**
## Parameters

For small amplitudes $Amp = 0.01$ at the input boundary condition the
system has analytical solution for $e1 < e2$, in this case the
spatial envelope of the flux rate behaves like $Q=Amp e^{-e2/2*x}$
[Wang et al., 2013].

## Initial conditions 

the initial conditions are $A=1$ and $Q=0$.
*/

event init (i = 0)
{
  foreach(){
    a[] = 1. ; 
    q[] = 0.0 ;
  }
  boundary({a,q});
}

int main() { 
  N = 512;
  e1 = 0.5 ;
  e2 = 0.1 ;
  omega = 1.;
  Amp = 0.01 ;
  DT = 1e-4;
  run();
}

/**
## Outputs

We print to the stderr the spatial profile of the flow rate $q$ and
the analytical solution. */
	
event printdata (t += 0.1; t <= 2.) {
  foreach()
    fprintf (stderr, "%g %g %.6f %g\n", x, q[], Amp*exp(-e2/2.*x), dt);
  fprintf (stderr, "\n\n");
}

/**
~~~gnuplot Comparison between numerical solution and linear theory
set yrange [0.008:]
set ylabel 'Q'
set xlabel 'x'
plot 'log' w l t 'numerical', 'log' u 1:3 w l t 'linear theory'
~~~

# Bibliography

*  [Wang et al., 2013](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/wang14.pdf)
*/
