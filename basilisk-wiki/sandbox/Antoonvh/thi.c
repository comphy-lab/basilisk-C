#define BGHOSTS 2
#include "grid/multigrid1D.h"
#include "higher-order.h"

#define GCA (sqrt(pi)/2.*(erf(x + (Delta/2.)) - erf(x - (Delta/2.)))/Delta)

scalar s[];

int main() {
  L0 = 15;
  X0 = Y0 = -L0/2;
  for (N = 16; N <= 1024; N *= 2) {
    init_grid (N);
    foreach()
      s[] = GCA;
    boundary({s});
    double err = 0;
    for (double xp = -2; xp <= 2; xp += 0.99)  
      err += fabs(exp(-sq(xp)) - interpolate_5(s, xp));
    printf ("%d %g\n", N, err);
  }
}

/**
~~~gnuplot L1 vs. N
set logscale x 2 
set logscale y
set xr [8:2048]
set grid
set size square
set xlabel 'N'
set ylabel 'L1'
plot 'out', 1e6*x**(-5)
~~~
 */
