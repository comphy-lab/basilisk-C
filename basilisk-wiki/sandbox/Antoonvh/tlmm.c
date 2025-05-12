/**
# Test of the BDF3 predictor corrector scheme

Using a single correction step. The test problem being"

$$s' = -5s,$$
$$s(0) = 1.$$


~~~gnuplot The 4-step solution seems unstable
set logscale x 2
set logscale y
set grid
set size square
set xlabel 'steps'
set ylabel 'error'
set xr[2:1024]
set yr[10**(-9):1]
plot 'out' t 'Data', 9*x**(-3), '3rd order'
~~~

~~~gnuplot Evolution of the coarse solution
reset
set xlabel 'time'
set ylabel 's'
set grid
plot 'evol' w l , exp(-5*x)
~~~
*/
#include "lmm.h"

scalar s[];

double k = -5; //Stiffness 
void expo (scalar * sl, double t, scalar * dsl) {
  scalar s, ds;
  foreach() {
    for (s, ds in sl, dsl)
      ds[] = k*s[];
  }
}

int main() {
  init_grid (1);
  for (DT = 0.25; DT > 0.001; DT /= 2) 
    run();
}

event init (t = 0) {
  foreach()
    s[] = 1;
}

event stepper (i++, last) {
  double DS = DT;
  dt = dtnext (DS);
  LM_BDF ({s}, dt, expo, i, 1);
  if (DT == 0.25) {
    static FILE * evol = fopen ("evol", "w");
    if (i == 0)
      fprintf (evol, "0 1\n");
    foreach()
      fprintf (evol, "%g %g\n", t + dt, s[]); 
  }
}

event stop (t = 1) {
  foreach()
    printf ("%d %g\n", i, fabs(s[] - exp(k*t)));
  return 1;
}
