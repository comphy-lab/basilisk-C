/**
# 5-point 4th-order accurate conservative refinement based on ordinary least squares

~~~gnuplot it works
set grid
set size square
set logscale x 2
set logscale y
set xr [10:1500]
set yr [10**-9:1]
plot 'out' t 'error method 1','' u 1:3 t 'error method 2', 999*x**-4 lw 3 t 'fourth order'
~~~

 */

#include "grid/bitree.h"
#include "utils.h"

void refine_4_ole (Point point, scalar s) {
  double a = //rare coefficients
    -489./4480. * s[-2] +
    1153./2240. * s[-1] +
      17./35.   * s[]   +
     383./2240. * s[1]  +
    -279./4480. *s[2];
  foreach_child()
    if (child.x == -1)
      s[] = a;
    else
      s[] = 2*coarse(s, 0) - a;
}

void refine_4_ole2 (Point point, scalar s) {
  double a = //rare coefficients
    -489./4480. * s[-2] +
    1153./2240. * s[-1] +
      17./35.   * s[]   +
     383./2240. * s[1]  +
    -279./4480. * s[2];
  double b = 
    -489./4480. * s[2]  +
    1153./2240. * s[1]  +
      17./35.   * s[]   +
     383./2240. * s[-1] +
    -279./4480. * s[-2];
  double c = s[] - (a + b)/2.;
  foreach_child()
    if (child.x == -1)
      s[] = a + c;
    else
      s[] = b + c;
}

int main() {
  L0 = 10;
  X0 = -L0/2.;
  
  for (N = 8; N <= 1024; N *= 2) {
    init_grid (N);
    scalar s[], s2[], w[], w2[];
    s.prolongation = refine_4_ole;
    s2.prolongation = refine_4_ole2;
    foreach()
      s2[] = s[] = exp(-sq(x));
    wavelet (s , w);
    wavelet (s2, w2);
    printf ("%d %g %g\n", N, normf(w).avg, normf(w2).avg);
  }
}
