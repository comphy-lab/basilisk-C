/**
## A 4th-order accurate compact finite volume Poisson solver 

![source term](compact_poisson/b.png)

![Solution](compact_poisson/a.png)

~~~gnuplot Convergence
set xr [20:600]
set logscale x 2
set logscale y
set grid
plot 'out' u 1:2 t 'L1', '' u 1:3 t 'Li', 1e3*x**-4 t 'fourth order'
~~~

~~~gnuplot Iteration convergence is poor on trees
reset
set xr [20:600]
set yr [25:35]
set grid
set logscale x 2
set xlabel 'N'
plot 'out' u 1:4 t 'MG Cycles'
~~~

## note

The iteration convergence is excellent of equidistant grids, but not on trees. 
*/

#include "cp.h"
#include "utils.h"

double func (double x, double y) {
  return ((4*(sq(x) + sq(y)) - 4)*exp(-sq(x) - sq(y)));
}

double func_sol (double x, double y) {
  return exp(-sq(x) - sq(y));
}

scalar a[], b[];
a[left] = dirichlet (0);

int main() {
  L0 = 10;
  X0 = -4.765; Y0 = -5.123;
  TOLERANCE = 1e-9;
  for (N = 32; N <= 512; N *= 2) {
    init_grid (N);
    int l = depth();
    refine (sq(x) + sq(y) < 1 && level <= l);
    a.prolongation = refine_5th;
    b.prolongation = refine_5th;
    foreach() {
      b[] = Gauss6 (x, y, Delta, func);
      a[]= 0;
    }
    mgstats mg = poisson_compact_mg (a, b);
    double e = 0, em = -1;
    foreach(reduction (+:e) reduction (max:em)) {
      double el = fabs(Gauss6(x, y, Delta, func_sol) - a[]);
      e += el*sq(Delta);
      if (el > em)
	em = el;
    }
    
    printf ("%d %g %g %d %d\n",
	    N, e, em, mg.i, mg.nrelax);
  }
  output_ppm (b, file = "b.png", n = 400);
  output_ppm (a, file = "a.png", n = 400, min = 0, max = 1);
}
