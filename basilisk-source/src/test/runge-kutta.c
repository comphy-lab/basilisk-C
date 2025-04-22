/**
# Convergence of the Runge--Kutta solvers

We solve numerically
$$
\frac{du}{dt} = tu
$$
with the initial condition $u(0)=1$. The solution is $u(t)=e^{t^2/2}$.
*/

#include "runge-kutta.h"

/**
This function returns the right-hand-side of the equation i.e. $tu$. */

static void du (scalar * ul, double t, scalar * kl)
{
  scalar u = ul[0], k = kl[0];
  foreach()
    k[] = t*u[];
}

int main()
{
  init_grid (1);

  for (int order = 1; order <= 4; order *= 2)
    for (double dt = 1e-2; dt <= 8e-2; dt *= 2) {
      scalar u[];
      foreach()
	u[] = 1.; // the initial condition
      double emax = 0.;
      for (t = 0; t <= 2.; t += dt) {
	foreach() {
	  double e = fabs (u[] - exp(t*t/2.));
	  if (e > emax)
	    emax = e;
	  printf ("%g %g %g\n", t, u[], u[] - exp(t*t/2.));
	}
	runge_kutta ({u}, t, dt, du, order);
      }
      printf ("\n");
      fprintf (stderr, "%g %g %d\n", dt, emax, order);
    }
}

/**
~~~gnuplot Error convergence for different orders
set xlabel 'dt'
set ylabel 'error'
set logscale
set key top left
set ytics format '%.0e'
plot "< grep '1$' log" pt 7 t '', 15.*x t '15 dt',  \
     "< grep '2$' log" pt 7 t '', 4.*x*x t '4 dt^2', \
     "< grep '4$' log" pt 7 t '', x**4/2. t 'dt^4/2'
~~~
*/
