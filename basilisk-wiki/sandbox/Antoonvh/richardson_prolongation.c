/**
# Richardson extrapolation prolongation

Using parent and grand parent cell data:

~~~gnuplot
set logscale x 2
set logscale y
set xr [20:700]
set yr [0.0001:0.5]
set xlabel 'N'
set ylabel 'Error'
set size square
set grid
plot 'out' t 'L_1', 100*x**(-2)
~~~
*/
#include "grid/multigrid1D.h"

static inline void refine_richardson_1 (Point point, scalar v) {
  double val = v[];
  double valp = val;
  int c = (child.x > 0 ? -1 : 1);
  if (level > 0)
    valp = coarse (v);
  foreach_child() 
    v[] = val + c*(child.x > 0 ? -1 : 1)*(val - valp)/2.;
}

int main() {
  L0 = 10;
  X0 = -L0/2.;
  for (N = 32; N <= 512; N *= 2) {
    init_grid (N);
    scalar s[], w[];
    s.prolongation = refine_richardson_1;
    //s.prolongation = refine_injection;
    foreach()
      s[] = exp(-sq(x));
    wavelet (s, w);
    double wt = 0;
    foreach (reduction (+:wt))
      wt += Delta*fabs(w[]);
    printf ("%ld %g\n", grid->tn, wt);
  }
}
