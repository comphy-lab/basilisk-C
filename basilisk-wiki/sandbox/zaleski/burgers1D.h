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
    u[] = sin(2.*pi*x);
}

event logfile (i++)
{
  stats s = statsf (u);
  fprintf (stderr, "%g %d %g %g %g %.8f\n", t, i, dt, s.min, s.max, s.sum);
}

event outputfile (t <= 0.5; t += 0.1)
{
  foreach()
    fprintf (stdout, "%g %g\n", x, u[]);
  fprintf (stderr, "\n");
}

/**
# Results

The sine wave steepens and forms a shock near $t_s = 1/(2\pi) \approx 0.159$.
After shock formation the central scheme will develop oscillations, as expected
for a non-dissipative scheme.

~~~gnuplot Evolution of the Burgers field.
plot 'out' u 1:2 w l t "u"
~~~
*/

/**
# To do

# Usage
  - [burgers1D.c]().
*/
   
