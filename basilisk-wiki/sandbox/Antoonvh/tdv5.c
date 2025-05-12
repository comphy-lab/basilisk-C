/**
# 5th-order accurate TVD advection. 

Error for advection of a Gaussian pulse that propaged a $1000$ times
its (Gaussian) width

~~~gnuplot Error in Smooth solution
set logscale x 2
set logscale y
set xr [24:560]
set grid
set size square
set ylabel 'L_1 error'
set xlabel 'N'
plot 'out' u 1:2, 1e9*x**(-5)
~~~

~~~gnuplot Total variance of advected noise
set ylabel 'Total Variance'
set key top left
plot 'out' u 1:3 t 't = 0', '' u 1:4 t 't = 1000'
~~~
*/

#define RKORDER (4)
#include "grid/bitree.h"
#include "lsrk.h"

void adv (scalar * sl, scalar * dsl) {
  boundary (sl);
  scalar s, ds;
  for (s, ds in sl, dsl) {
    scalar dsdt[];
    foreach() 
      dsdt[] = -(s[1] - s[0])/Delta;
    boundary ({ds});
    scalar rhs[];
    foreach() {
      ds[] = ((-s[-2] + 7.*(s[-1] + s[]) - s[1])/12.);
      rhs[] = 251./152.*dsdt[-1] + 1./8.*dsdt[-2];
    }
    for (int it = 0; it < 8; it++) {
      foreach() {
	double das = rhs[] - 413./456.*ds[-1] + 23./152.*ds[1] - 5./228*ds[2];
	ds[] = (ds[] + 8.*das)/9.; // 8/9th under relaxation
      }
      boundary ({ds});
    }
  }
}

double fun (double x) {
  return (exp(-sq(x)));
}
  
scalar s[], n[]; // Smooth and noisy fields

int main() {
  periodic (left);
  L0 = 20;
  X0 = -L0/2;
  double CFL = 0.36;
  for (N = 32; N <= 512; N *= 2) {
    DT = CFL*L0/(N); 
    run();
  }
}

double V1;
event init (t = 0) {
  foreach() {
    s[] = fun(x - Delta/2.);
    n[] = noise();
  }
  boundary ({n});
  foreach()
    V1 += fabs(n[1] - n[0]);
}

event advance (i++, last) {
  dt = dtnext (DT);
  A_Time_Step ({s, n}, dt, adv);
}

event end (t = 1000.0) {
  double V = 0, e = 0;
  foreach() {
    V += fabs(n[1] - n[0]);
    e += Delta * fabs(s[] - fun(x - Delta/2));
  }
  printf ("%d %g %g %g\n", N, e, V1, V);
  return 1;
}  
