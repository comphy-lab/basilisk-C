/**					
# Compact upwind advection

The compact upwind advection scheme of Van Hooft and Popinet (2022) is
tested for a smooth non-analytic function and a Gaussian.
*/
double smooth_bump (double x) {
  if (fabs(x) >= 1)
    return 0;
  return exp(-1./(1 - sq(x)) + 1);
}
/**
## Visual inspection:

~~~gnuplot Initial and final solutions for N = 256
set xlabel 'x'
set xr [-4:3]
set yr [ -.1:1.1]
set grid
plot 'log' w l t 'Smooth t = 0', '' u 1:3 w l t 'Smooth t = 1',  '' u 1:4 w l t 'Gaussian t = 0', '' u 1:5 w l t 'Gaussian t = 1' 
~~~

~~~gnuplot You may inspect the Smooth function near the transition at x= 1
set xr [0.5:1.1]
set yr [-0.1:0.6]
set grid
plot (x >= 1) ? 0 : exp(-1/(1 - x*x) + 1) w l lw 2
~~~

## Convergence

~~~gnuplot Convergence for the Gaussian function is at a fourth-order rate
set logscale x 2
set logscale y
set grid
set size square
set xr [16:2048]
set yr [1e-8:0.1]
plot 'out' u 1:4  t 'L2 Gaussian', '' u 1:5 t 'L_{inf} Gaussian', 1e5*x**(-4) t 'Fourth order'
~~~

~~~gnuplot Convergence for the smooth function is only eventually at a fourth-th order rate.
set logscale x 2
set logscale y
set grid
set size square
set xr [16:17000]
set yr [1e-8:1]
plot 'out'  t 'L2 Smooth', '' u 1:3 t 'L_{inf} smooth', 1e8*x**(-4) t 'Fourth order'
~~~

Thank you `smooth_bump` for having such sharp features, giving rise to this critical note.
*/
#include "grid/multigrid1D.h"
#define RKORDER 4
#include "higher-order.h"

#include "lsrk.h"

vector u[];

void advection (scalar * sl, scalar * dsl) {
  scalar s, ds;
  for (s, ds in sl, dsl) {
    vector grad[];
    compact_upwind ({s}, {grad}, u);
    foreach()
      ds[] = -u.x[]*grad.x[];
  }
}

scalar smooth[], Gauss[];
scalar * tracers = {smooth, Gauss};

int main() {
  periodic (left);
  L0 = 10;
  X0 = -6;
  for (N = 32; N <= 10000; N *= 2) {
    DT = L0/N;
    run();
  }
}

event init (t = 0) {
  dt = dtnext (DT);
  foreach() {
    u.x[] = 1;
    smooth[] = smooth_bump(x + 1);
    Gauss[] = exp(-sq(x + 1));
  }
}

event advance (i++, last) {
  dt = dtnext (DT);
  A_Time_Step (tracers, dt, advection);
  boundary (tracers);
}

event stop (t = 1) {
  double L2s = 0, Lis = 0;
  double L2g = 0, Lig = 0;
  
  foreach(reduction(+:L2s) reduction(max:Lis) reduction(+:L2g) reduction(max:Lig)) {
    if (N == 256)
      fprintf (stderr, "%g %g %g %g %g\n", x, smooth_bump(x + 1), smooth[], exp(-sq(x + 1)), Gauss[]);
    L2s += Delta*sq(smooth[] - smooth_bump(x));
    L2g += Delta*sq(Gauss[] - exp(-sq(x)));
    if (fabs(smooth[] - smooth_bump(x)) > Lis)
	Lis = fabs(smooth[] - smooth_bump(x));
    if (fabs(Gauss[] - exp(-sq(x))) > Lig)
      Lig = fabs(Gauss[] - exp(-sq(x)));

  }
  printf ("%d %g %g %g %g\n", N, sqrt(L2s), Lis, sqrt(L2g), Lig);
}

/**
## Reference

van Hooft, J. Antoon, and Stéphane Popinet. "A fourth-order accurate adaptive solver for incompressible flow problems." Journal of Computational Physics 462 (2022): 111251.
 */