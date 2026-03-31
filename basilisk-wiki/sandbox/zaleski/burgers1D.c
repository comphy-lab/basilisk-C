/**
# Burgers equation — sine wave initial condition, explicit central scheme in 1D.

Solves the inviscid Burgers equation
$$
\partial_t u + \partial_x\!\left(\frac{u^2}{2}\right) = 0
$$
with a sine-wave initial condition $u(x,0) = \sin(2\pi x)$, which produces
a shock at time $t_s = 1/(2\pi) \approx 0.159$.

See [centralscheme.h]() for the time integration scheme.
*/

#include "grid/cartesian1D.h"
#include "centralscheme.h"

#define LEVEL 8

int main()
{
  origin (-0.5, 0.);
  init_grid (1 << LEVEL);
  DT = 1e-3;
  run();
}

event init (i = 0)
{
  foreach()
    u[] = sin(2.*pi*x/L0);
}

event logfile (i++)
{
  stats s = statsf (u);
  fprintf (stderr, "%g %d %g %g %g %.8f\n", t, i, dt, s.min, s.max, s.sum);
}

event outputfile (t <= 0.2; t += 0.05)
{
  foreach()
    fprintf (stdout, "%g %g\n", x, u[]);
  fprintf (stderr, "\n");
  fprintf (stdout, "\n\n");
}


event energy (i++)
{
  double ke = 0.;
  foreach()
    ke += 0.5*sq(u[])*Delta;
  static FILE * fpe = NULL;
  if (fpe == NULL)
    fpe = fopen ("energy", "w");
  fprintf (fpe, "%g %g\n", t, ke);
}

/**
# Results

The sine wave steepens and forms a shock near $t_s = 1/(2\pi) \approx 0.159$.
**Even before shock formation** the central scheme develops oscillations (as expected for a non-dissipative scheme?).

~~~gnuplot Evolution of the Burgers field.
set key top left
plot for [i=0:4] 'out' index i u 1:2 w l t sprintf("t = %3.1f", i*0.05)
~~~

The kinetic energy $E_K = \int_0^{L_0} \frac{u^2}2 dx$ is plotted. The energy is slowly increasing at first (perhaps an effect of the time integration) then decreasing once the shock forms (perhaps because of energy dissipation inside the shock). 

~~~gnuplot Kinetic energy as a function of time. 

plot 'energy' u 1:2 w l t "E_k"
~~~
*/
