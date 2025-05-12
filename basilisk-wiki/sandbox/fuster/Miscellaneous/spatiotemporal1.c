/**                                                                                                                                                                                         
# 1D equivalent of spatio-temporal problem
*/

#include "grid/multigrid1D.h"
#include "run.h"
#include "diffusion.h"

scalar f[];

f[left] = dirichlet(1);
f[right] = dirichlet(2);

FILE * fptrans;

int main() {

  X0 = 0;


  fptrans = fopen("error.dat", "w");

  for (int N =32; N <= 1024; N *= 2) {
      init_grid (N);
      run();
      free_grid();
  }

  fclose(fptrans);
}

event output (t += 0.1; t <= 1) {

  foreach()
    printf ("%g %g %g \n", x, f[], t);
  printf ("\n");

}


event integration (i++) {

  dt = dtnext (1./N);
  scalar r[];  
  foreach()
    r[] = 2.; 
  diffusion (f, dt, r = r); 

}

  
double errtrans = 0.;
double errtransL1 = 0.;

event transient ( i++ ) {

  if ( t > 0.01 - 1.5/N && t < 0.01 + 1.5/N && t > 0.) {

    double tmp = 0.;
    double tmpL1 = 0.;
    foreach () {
      if (x < 0.1) {
        tmp = max(tmp, fabs(erf(x/(2.*sqrt(t))) -1. + f[]));
        tmpL1 += fabs(erf(x/(2.*sqrt(t))) -1. + f[])*Delta;
      }
      if (y > 0.9) {
        tmp = max(tmp, fabs(erf((1.-x)/(2.*sqrt(t))) - (2. - f[])/2.));
        tmpL1 += fabs(erf(x/(2.*sqrt(t))) -1. + f[])*Delta;
      }
    }
    fprintf(fptrans, "%i %g %g \n", N, tmp, tmpL1);

    errtrans = max(errtrans, tmp);
    errtransL1 = max(errtransL1, tmpL1);
  }

}

event error ( t = 1) {

  double err = 0., errmax = 0.;
  foreach () {
    err   += fabs( f[] - ( 1. + 2*x - sq(x) ) )*Delta;
    errmax = max(errmax, fabs( f[] - ( 1. + 2*x - sq(x) ) ));
  }
  fprintf(stderr, "%i %g %g %g %g %g \n", N, perf.t, err, errmax, errtrans, errtransL1);

}

/**
~~~gnuplot
plot 'out' u 1:2:3 w l palette, 1 - x**2 + 2.*x t 'EXACT' w l lw 3 lc 0
~~~

~~~gnuplot
set log xy
set key below
set xlabel 'N'
plot 'log' u 1:4 t 'steady L_{max}' w lp, 'error.dat' u 1:2 t 'transient L_{max}' w p, 'log' u 1:3 t 'steady L1' w p, 'error.dat' u 1:3 t 'transient L1' w p
~~~

*/
