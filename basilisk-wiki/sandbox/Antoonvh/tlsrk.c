/**
# A test for the Low-Storage RK integrator

$$s_t = s,$$
$$s(0) = 1.$$


~~~gnuplot
set xr [0.5:16000]
set yr [1e-10:1e5]
set logscale xy
set xlabel 'Error'
set ylabel 'Steps'
set size square
set grid
plot 'out' t 'data', 1e6*x**(-4) t '4th order'
~~~
*/
#define RKORDER 4
#include "lsrk.h"

void df_is_f (scalar * al, scalar * bl) {
  scalar a, b;
  foreach() {
    for (a, b in al, bl) 
      b[] =  a[];
  }
}

double tend = 10;
scalar s[];

int main () {
  N = 1;
  for (DT = tend; DT > 0.001; DT /= 2)
    run();
}

event init (t = 0) {
  foreach()
    s[] = 1;
}

event step (i++, last) {
  dt = dtnext (DT);
  A_Time_Step ({s}, DT, df_is_f);
}
  
event stop (t = tend) {
  foreach()
    printf ("%d %g\n", (int) (tend/DT + 0.5), fabs(s[] - exp(tend)));
  return 1;
}


