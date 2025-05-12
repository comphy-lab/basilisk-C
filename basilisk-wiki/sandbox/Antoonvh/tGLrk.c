/**
# A test for the 6-th order implicit RK integrator

$$s_t = s,$$
$$s(0) = 1.$$

~~~gnuplot
set xr [0.5:16000]
set yr [1e-12:1e5]
set logscale xy
set xlabel 'Steps'
set ylabel 'Error'
set size square
set grid
plot 'out' t 'data', 1e7*x**(-6) t '6th order'
~~~
*/
#define RKORDER 6
#include "GLrk.h"
double Tol = 1e-4;

void df_is_f (scalar * al, scalar * bl) {
  foreach() {
    scalar a, b;
    for (a, b in al, bl)
      b[] = a[];
  }
}

double tend = 10;
scalar s[];

int main () {
  N = 1;
  for (DT = 5; DT > 0.001; DT /= 2)
    run();
}

event init (t = 0) {
  foreach()
    s[] = 1;
}

event step (i++, last) {
  dt = dtnext (DT);
  A_Time_Step ({s}, DT, df_is_f, Tol);
}

event stop (t = tend) {
  double e = 1;
  foreach()
    e = fabs(s[] - exp(tend));
  printf ("%d %g\n", (int) (tend/DT + 0.5), e);
  Tol = e/(pow(2,RKORDER)*100*tend);
  return 1;
}
