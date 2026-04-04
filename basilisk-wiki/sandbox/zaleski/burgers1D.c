/**
# Burgers equation — sine wave initial condition, explicit central schemes in 1D.

Solves the inviscid Burgers equation
$$
\partial_t u + \partial_x\!\left(\frac{u^2}{2}\right) = 0
$$
with a sine-wave initial condition $u(x,0) = \sin(2\pi x)$, which produces
a shock at time $t_s = 1/(2\pi) \approx 0.159$.

See [centralscheme.h](), [arakawa.h]() and [midpoint-skew-symmetric.h]()
for the time integration schemes.

This version ([midpoint-skew-symmetric.h]() ) uses a midpoint skew symmetric scheme. The midpoint is obtained by iterations, in this case four iterations are used.  With four iterations as below the energy has an escursion $E_k(t) - E_k(0)$ of $10^{-7}$. 
With one iteration the energy has an excursion of $4 10^{-5}$. 
With the explicit skew-symmetric scheme  [arakawa.h]() the energy has an excursion of $10^{-4}$. 

*/

#include "grid/cartesian1D.h"
// #include "centralscheme.h"

#define NITER 4
#include "midpoint-skew-symmetric.h"

//#include "arakawa.h"

#define LEVEL 8

int main()
{
  origin (-0.5);
  // periodic (right);
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

event outputfile (t <= 0.5; t += 0.05)
{
  foreach()
    fprintf (stdout, "%g %g\n", x, u[]);
  fprintf (stderr, "\n");
  fprintf (stdout, "\n\n");
}

event energy (i++)
{
  double ke = -0.25;
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

~~~gnuplot Evolution of the Burgers field after shock formation.
set key top left
plot for [i=5:9:2] 'out' index i u 1:2 w l t sprintf("t = %3.1f", i*0.05)
~~~

The **kinetic energy** $E_K = \int_0^{L_0} \frac{u^2}2 dx$ is plotted.

~~~gnuplot Kinetic energy as a function of time.
set grid
plot 'energy' u 1:2 w l t "E_k(t) - E_k(0)"
~~~

~~~gnuplot Kinetic energy at early times.
set xrange [0:0.05]
plot 'energy' u 1:2 w l t "E_k"
~~~

~~~gnuplot Evolution of the residual
set key top left
plot for [i=0:500:200] 'residual' index i u 1:2 w l t sprintf("t = %3.1f", i*0.05)
~~~


*/
