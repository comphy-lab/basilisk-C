/**
# Viscous decay of a Sinusoidal flow profile

$$u_y = \mathrm{sin}(x)e^{-\nu t}$$

![The test](decay/uy.mp4)

~~~gnuplot Convergence
set xr[4:256]
set grid
set logscale x 2
set logscale y
set xlabel 'N'
set ylabel 'L_1'
set size square
plot 'out' t 'data', 1e3*x**(-4)
~~~
 */
#include "nsf4t.h"
scalar * tracers = NULL;
double muv = 0.1;

double uy (double x, double y) {
  return (sin(x)*exp(-muv*t));
}
    
int main() {
  foreach_dimension()
    periodic (left);
  L0 = 2*pi;
  const scalar muc[] = 0.1;
  nu = muc;
  for (N = 8; N <= 128; N *= 2)
    run();
}

event init (t = 0) {
  CFL = 0.2; // Stability at N = 128
  foreach_face(y)
    u.y[] = Gauss6_y (x, y, Delta, uy);
}

event mov (i += 5) {
  if (N == 128)
    output_ppm (u.y, file = "uy.mp4", n = 300, min = -1.1, max = 1.1);
}

event stop (t = 10) {
   double e = 0;
   foreach_face(y)
     e += dv() * fabs(u.y[] - Gauss6_y(x, y, Delta, uy));
   printf ("%d %g\n", N, e);
}
