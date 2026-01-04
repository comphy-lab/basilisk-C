/**
# Compact Poisson-Helmholtz problem.

We combine the first-order implicit time discretization of
`diffusion.h` with *the* 5-point implicit eight-order second
derivative (compact scheme). A "double implicit" approach.

~~~gnuplot Diffusion $\Delta t = 2$ and $\kappa = 0.5$
set xlabel 'x'
set ylabel 's'
plot 'out' t 'After a time step', '' u 1:3 t 'Before'
~~~

Using the soltution at $t = t_0 + \Delta t$ we apply a "forward", step with
negative $\Delta t$, using the standard coefficients of eight-order accurate
second derivative to check if we arrive at the initial solution.

~~~gnuplot Error is within a set tolerance
set yr [0:2e-9]
set xlabel 'x'
set ylabel 'Error'
plot 'log' u 1:3 t 'error', 1e-9 t 'Tolerance'
~~~
 */
#include "grid/multigrid1D.h"
#include "solve.h"

//1D Coefs as in Lele (1992)
double alpha = 344./1179, beta, a, b, c;

void compute_constants_1D() {
  beta = (38*alpha - 9)/214.;
  a = (696 - 1191*alpha)/428.;
  b = (2454*alpha - 294)/535.;
}

int main() {
  compute_constants_1D();
  periodic (left);
  TOLERANCE = 1e-9;
  L0 = 20;
  X0 = -L0/2.;
  init_grid (128);
  scalar s[], sn[];
  foreach() {
    s[] = exp(-sq(x));
    sn[] = 0;
  }
  // Implicit time step
  double dt = 2;
  double kappa = 0.5;
  solve (sn, kappa*((a*(sn[-1] - 2*sn[] + sn[1]) + b*(sn[-2] - 2*sn[] + sn[2])/4.)/sq(Delta)) +
	 -1./dt*(beta*(sn[2] + sn[-2]) + alpha*(sn[1] + sn[-1]) + sn[]) ,
	 -1./dt*(beta*(s[2] + s[-2]) + alpha*(s[1] + s[-1]) + s[]));
  foreach() {
    printf ("%g %g %g\n", x, sn[], s[]);
    s[] = sn[];
    sn[] = 0.;
  }
  // Compute second derivative
  solve (sn, beta*(sn[2] + sn[-2]) + alpha*(sn[1] + sn[-1]) + sn[] ,
	 (a*(s[-1] - 2*s[] + s[1]) + b*(s[-2] - 2*s[] + s[2])/4.)/sq(Delta));
  // Andvance backward in time to verify that we have used the double derivative
  foreach()
    s[] -= dt*kappa*sn[];
  foreach()
    fprintf (stderr, "%g %g %g\n", x, s[], fabs(s[] - exp(-sq(x))));
}