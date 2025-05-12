/**
# The steady 3D flow of Antuono

The test case is too simple as we get 5th-order convergence with our
4th-order accurate formulations.

~~~gnuplot Error convergence
set xr [4:128]
set logscale x 2
set logscale y
set xlabel 'N'
set ylabel 'Error'
plot 'out' u 1:2 t 'L_1', '' u 1:3 t 'Max', 1e4*x**(-5)
~~~

See also [here](antuono.c)
 */
#include "grid/octree.h"
#include "nsf4t.h"
scalar * tracers = NULL;

foreach_dimension()
double Antuono_x (double x, double y, double z) {
  double a = 4*sqrt(2)/(3*sqrt(3));
  coord c = {x, y ,z};
  return a*(sin(c.x - 5.*pi/6)*cos(c.y - pi/6)*sin(c.z) -
	    cos(c.z - 5.*pi/6)*sin(c.x - pi/6)*sin(c.y));
}

int main() {
  foreach_dimension()
    periodic (left);
  L0 = 2*pi;
  for (N = 8; N <= 64; N *= 2)
    run();
}

event init (t = 0) {
  TOLERANCE = 1e-5;
  foreach_face()
    u.x[] = Gauss6_x(x, y, z, Delta, Antuono_x);
}

event stop (t = 1) {
  double e = 0, em = -1;
  foreach_face() {
    double el = fabs(u.x[] - Gauss6_x(x, y, z, Delta, Antuono_x));
    e += dv()*el;
    if (el > em)
      em = el;
  }
  printf ("%d %g %g\n", N, e, em);
  return 1;
}

