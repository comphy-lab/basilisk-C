/**
# Using `solve.h` to evaluate a compact scheme

We try the 4th-order Pade scheme:

~~~gnuplot Convergence is good, upto the TOLERANCE
set xr [8 : 2048]
set logscale x 2
set logscale y
set grid
set xlabel 'N'
set ylabel 'error'
plot 'out' u 1:2 t 'L1', '' u 1:3 t 'L_i', 2e3*x**-4 t 'Fourth order', '' u 1:6 t 'Residual'
~~~

The convergence rate is quite good

~~~gnuplot Solver statistics
reset
set xr [8 : 2048]
set yr [0:10]
set xlabel 'N'
set logscale x 2
set grid
plot 'out' u 1:4 t 'MG cycles', 'out' u 1:5 t 'relaxations'
~~~
 */
#include "solve.h"

int main() {
  TOLERANCE = 1e-6; // Much larger than the minimum error
  L0 = 10;
  X0 = Y0 = -L0/2;
  for (N = 16; N <= 1024; N *= 2) {
    init_grid (N);
    scalar s[], ds[], rhs[];
    foreach() {
      s[] = exp(-sq(x) - sq(y));
      ds[] = 0;
    }
    /**
       Pade scheme: 
       
       $$\frac{1}{4}f'_{i-1} + f'_{i} + \frac{1}{4}f'_{i+1} =
       \frac{3}{2} \frac{f_{i+1}-f_{i-1}}{2\Delta}.$$
    */
    foreach() 
      rhs[] = 3.*(s[1] - s[-1])/(4.*Delta);
    mgstats stats = solve (ds, ds[-1]/4. + ds[] + ds[1]/4., rhs[]); 
    double e = 0, em = -1;
    foreach(reduction (+:e) reduction (max:em)) {
      double el = fabs(-2*x*exp(-sq(x) - sq(y)) - ds[]);
      if (el > em)
	em = el;
      e += sq(Delta)*el;
    }
    printf ("%d %g %g %d %d %g\n", N, e, em, stats.i, stats.nrelax, stats.resa);
  }
}
