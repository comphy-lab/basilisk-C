/**
# Self-similar pinch-off

Here we check whether we can recover the self-similar scalings with
time of the minimum radius and maximum axial velocity. */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

double R = 0.2, epsilon = 0.1;

int main() {
  rho2 = 0.01;
  f.sigma = 1.;
  N = 256;
  run();
}

event init (t = 0) {
  fraction (f, R*(1. + epsilon*cos(pi*x)) - y);
}

event adapt (i++) {
  adapt_wavelet ({f}, (double[]){0}, 8); // t < 0.7 ? 8 : 9);
}



event logfile (i++; t<=1) {
  scalar hy[];
  foreach() {
    if (f[] > 1e-3 && f[] < 1 - 1e-3) {
      coord m = mycs (point, f), fc;
      double alpha = plane_alpha (f[], m);
      plane_area_center (m, alpha, &fc);
      hy[] = y + Delta*fc.y;
    }
    else
      hy[] = nodata;
  }
  stats s = statsf (hy);
  fprintf (stderr, "%g %g %g %g\n", t, s.min, s.max, statsf(u.x).max);
}

/**
~~~gnuplot Evolution of the minimum radius with time
set xlabel 'Time'
set ylabel 'Rmin'
set logscale
set grid
t0=0.7541
plot [0.001:]'log' u (t0-$1):2 w l, x**(2./3.)/1.7
~~~

~~~gnuplot Evolution of the maximum velocity with time
set ylabel 'Maximum velocity'
plot [0.001:]'log' u (t0-$1):4 w l, x**(-1/3.)/7.
~~~

## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/plateau.html)
*/
