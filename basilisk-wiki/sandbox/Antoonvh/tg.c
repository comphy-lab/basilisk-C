/**
# Taylor-Green vortices

![Advecting vortices](tg/o.mp4)

~~~gnuplot L1 Error
set xr [8: 256]
set xlabel 'Spatial resolution (N)'
set ylabel 'Error'
set size square
set grid
set logscale xy 2
plot 'out' u 1:2, 1e3*x**(-4) t '4th order'
~~~
 */
#define RKORDER (4)
#define VERT_4 1
#define FACE_4 1
#define ADV_4 1
#define PROJECT_4 1
#include "nsf2.h"

int main() {
  origin (-0.5, -0.5);
  foreach_dimension()
    periodic (right);
  for (N = 16; N <= 128; N *= 2)
    run();
}

double u_x (double x, double y) {
  return -cos(2.*pi*x)*sin(2.*pi*y) + 1.;
}

double u_y (double x, double y) {
  return  sin(2.*pi*x)*cos(2.*pi*y) + 0.5;
}

event init (t = 0) {
  foreach_face() 
    u.x[] = Gauss6_x (x, y, Delta, u_x);
}

event mov (i += 5) {
  scalar omg[];
  vorticityf (u, omg);
  output_ppm (omg, file = "o.mp4", n = 300);
}

event error (t = 2) {
  double e = 0;
  foreach_face() 
    e += fabs(Gauss6_x (x, y, Delta, u_x) - u.x[]);
  printf ("%d %g\n", N, e/sq(N));
}
