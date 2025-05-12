/**
# A fluid over an instantanious moving plate. 

We test the accuracy of the centered solver for a non-steady *Stokes* flow.
*/
#include "navier-stokes/centered.h"

#define sol(y) (1 - erf(y/(2.*pow(t, 0.5))))
u.t[bottom] = dirichlet(1.);
u.t[top] = dirichlet(t > 0 ? sol(y) : 0);

int main(){
  L0 = 10;
  periodic(left);
  const face vector muc[] = {1., 1.};
  mu = muc;
  TOLERANCE = 10E-5;
  DT = 0.001; // Errors due to time integration should remain very small
  for (N = 8; N <= 64; N *= 2){
    init_grid(N);
    run();
  }
}

event init (t = 0){
  foreach()
    u.x[] = 0;
}

event stop(t = 1){
  static FILE * fp = fopen("convergence", "w");
  double e = 0;
  foreach(){
    e += fabs(u.x[] - sol(y)) * sq(Delta);
  }
  fprintf(fp, "%d\t%g\n", N, e);
  fflush(fp);
  return(1);
}
/**
~~~gnuplot
set xr [5 : 127]
set logscale xy 2
set xlabel 'N'
set ylabel 'L1 Error'
set size square
set key box
plot 'convergence' u 1:2 t 'Data' , 10*x**(-2) w l t 'second order'
~~~
*/



