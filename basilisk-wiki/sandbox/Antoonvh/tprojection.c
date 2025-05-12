/**
## 4th order projection on a tree grid

Following the [Hodge decomposition example](../jmf/Hodge/hodge.c) of
J.M. Fullana, we test the accuracy of the 4th order projection method
on trees.

![Refine tree level](tprojection/lev.mp4)

~~~gnuplot L1 error
set xr [5 : 400]
set yr [5e-6:2]
set grid
set logscale x 2
set logscale y
set size square
set xlabel 'N_{eff}'
set ylabel 'L_1'
plot 'out' u 2:3 t 'data', 1e5*x**(-4)
~~~

~~~gnuplot Maximum divergence
reset
set xr [8 : 1024]
set yr [5e-12:1.1e-10]
set grid
set logscale x 2
set size square
set xlabel 'N'
set ylabel 'Max divergence'
plot 'out' u 1:6 t 'data', 1e-10 w l lw 2 t 'Tolerance'
~~~
 */

#include "nsf4t.h"
scalar * tracers = NULL;

#define GAUSSIAN (exp(-sq((x)) - sq((y))))

double us_x (double x, double y) {
  return (-2*y*GAUSSIAN);
}

double us_y (double x, double y) {
  return ( 2*x*GAUSSIAN);
}

double ui_x (double x, double y) {
  return ( 2*x*GAUSSIAN);
}

double ui_y (double x, double y) {
  return ( 2*y*GAUSSIAN);
}

int main() {
  foreach_dimension()
    periodic (left);
  p.prolongation = refine_4th;
  L0 = 20;
  X0 = -9.9; Y0 = -10.05;
  TOLERANCE = 1e-10;
  for (int l = 4; l <= 9; l++) {
    init_grid (1 << l);
    refine (sq(x) + sq(y) < sq(1) && level <= l);
    unrefine (sq(x) + sq(y - 0.5) > sq(2) && level >= l - 1);
    foreach_face() 
      u.x[] = Gauss6_x(x, y, Delta, us_x) + Gauss6_x(x, y, Delta, ui_x);
    boundary_flux ({u});
    mgstats mg = project (u, p);
    double e = 0;
    foreach_face() 
      e += dv()*fabs(u.x[] - Gauss6_x(x,y,Delta, us_x));
    scalar div[], lev[];
    foreach() {
      lev[] = level;
      div[] = 0;
      foreach_dimension()
	div[] += u.x[1] - u.x[];
      div[] /= Delta;
    }
    output_ppm (lev, file = "lev.mp4", n = 300, min = 0, max = 10.1, opt = "-r 2");
    stats ds = statsf (div);
    printf ("%d %g %g %d %d %g %g\n", N, sqrt(grid->tn), e, mg.i, mg.nrelax, ds.max, mg.resa);
  }
}
  
  
