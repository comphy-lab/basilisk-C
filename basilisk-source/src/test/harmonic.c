/**
# Very simple harmonic analysis */

#include "run.h"
#include "harmonic.h"

scalar s[], e[];

int main()
{
  DT = HUGE;

  /**
  We use a single grid point. */
  
  init_grid (1);
  run();
}

/**
This is the signal we want to decompose into harmonics. We add a 10%
noise to make things more interesting. */

event timestep (i++)
{
  foreach()
    s[] = 1. + sin(t) + 2.*cos(2.*t) + 0.1*noise();
}

/**
The harmonic decomposition itself. We store the error in `e`. The
inputs are the harmonic frequencies 1 and 2. Be careful to terminate
the array with zero. */

event harmonic (i++)
{
  harmonic_decomposition (s, t, {1, 2, 0}, e);
}

/**
We check whether there has been enough data to perform the inversion,
and then print each of the decomposition coefficients i.e.
$$
Z + \sum_i A_i \cos(\omega_i t) + B_i \sin(\omega_i t)
$$
*/

event output (t += 0.5; t <= 2.*pi)
{
  if (s.harmonic.invertible) {
    scalar Z = s.harmonic.Z, A0 = s.harmonic.A[0], B0 = s.harmonic.B[0],
      A1 = s.harmonic.A[1], B1 = s.harmonic.B[1];
    foreach()
      fprintf (stderr, "%g %g %g %g %g %g %g %g\n", t, s[], Z[], A0[], B0[], A1[], B1[], e[]);
  }
}

/**
As expected, the fit improves with time, as the random noise is
averaged out by the least-mean-square fitting of the decomposition.

~~~gnuplot Convergence of the coefficients of the harmonic decomposition
set xlabel 'Time'
set grid
Z = real(system("tail -n1 log | cut -d' ' -f3"))
A0 = real(system("tail -n1 log | cut -d' ' -f4"))
B0 = real(system("tail -n1 log | cut -d' ' -f5"))
A1 = real(system("tail -n1 log | cut -d' ' -f6"))
B1 = real(system("tail -n1 log | cut -d' ' -f7"))
set key bottom left
plot 'log' u 1:2 t 'data', '' u 1:3 t 'Z', \
    '' u 1:4 t 'A0', '' u 1:5 t 'B0', \
    '' u 1:6 t 'A1' lt 12, '' u 1:7 t 'B1', \
    '' u 1:8 t 'e', \
    Z + A0*cos(x) + B0*sin(x) + A1*cos(2.*x) + B1*sin(2.*x) t 'harmonic fit'
~~~
*/
