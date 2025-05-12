/**
# A Viscous boundary layer

It is resolved with 4th order accuracy.

![`u.x` for `N = 64`](visc_boun/ux.mp4)

~~~gnuplot
set xr [4:256]
set size square
set grid
set xlabel 'N'
set ylabel 'error'
set logscale x 2
set logscale y
plot 'out' t 'L_1', '' u 1:3 t 'Max', 1e3*x**(-4)
~~~
 */
#define NOSLIP_TOP (1)
#include "nsf4t.h"
scalar * tracers = NULL;

double t0 = 1;
double sol (double x, double y) {
  return (erf((Y0 + L0 - y)/(2.*sqrt(t + t0))));
}

int main() {
  L0 = 15;
  periodic (left);
  const scalar mu[] = 1.;
  DT = 0.01;
  nu = mu;
  for (N = 8; N <= 128; N *= 2)
    run();
}

event init (t = 0) {
  DI = .1;
  TOLERANCE = 1e-10;
  foreach_face(x)
    u.x[] = Gauss6_x(x, y, Delta, sol);
}

event mov (i += 2) {
  if (N == 64)
    output_ppm (u.x, file = "ux.mp4", n = 300, min = -0.1, max = 1.1);
}

event stop (t = 1) {
  double e = 0, em = -1;
  foreach_face(x) {
    double el = fabs(u.x[] - Gauss6_x(x, y, Delta, sol));
    e += sq(Delta)*el;
    if (el > em)
      em = el;
  }
  printf ("%d %g %g\n", N, e, em);
}
