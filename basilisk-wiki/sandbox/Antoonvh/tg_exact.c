/**
# Taylor-Green vortices

The solver appears exact when some deliberate mistakes are made.

![Advecting vortices](tg_exact/o.mp4)

~~~gnuplot L1 Error (mind the axis)
set xr [8: 256]
set yr [1e-16: 1e-13]
set xlabel 'Spatial resolution (N)'
set ylabel 'Error'
set size square
set grid
set logscale x 2
plot 'out' u 1:2
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
  return -cos(2.*pi*x)*sin(2.*pi*y);
}

double u_y (double x, double y) {
  return  sin(2.*pi*x)*cos(2.*pi*y);
}

event init (t = 0) {
  foreach_face() 
    u.x[] = u_x (x, y);
}

event mov (i += 5) {
  scalar omg[];
  vorticityf (u, omg);
  output_ppm (omg, file = "o.mp4", n = 300);
}

event error (t = 2) {
  double e = 0;
  foreach_face() 
    e += fabs(u_x (x, y) - u.x[]);
  printf ("%d %g\n", N, e/sq(N));
}
