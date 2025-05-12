/**
# Convergence for the first derivative of a function

We start with the previous code 
to compute the *convergence rate* of the approximation. */                                                                                            

#include "grid/cartesian1D.h"
#include "utils.h"

int main() {

  double sigma = 0.1;
  double x0 = 0.3;

  /** We replace the constant value of *n* with a loop going through
  10, 20, 40, 80 and 160. */
  for (int n = 10; n <= 160; n *= 2) {
    init_grid (n);

    scalar f[], dfexact[];
    foreach() {
      f[] = erf((x-x0)/sigma);
      dfexact[] = 2./sigma/sqrt(M_PI)*exp(-pow((x-Delta/2.-x0)/sigma,2));
    }   
    
    scalar df[];
    foreach()
      df[] = (f[] - f[-1])/Delta;

    scalar e[];
    foreach()
      e[] = fabs(df[] - dfexact[]);
    
    printf ("%d %g\n", n, normf(e).max);
  }
}

/** After the program has finished we plot the error as a function of *1/Dx*

~~~gnuplot
set output 'convergence.png'
set logscale
set xlabel '1/dx'
set ylabel 'normf(e).max'
fit a*x+b 'out' u (log($1)):(log($2)) via a,b
plot [5:320]'out' pt 7 t '', exp(b)*x**a t sprintf("%.0f/n^{%4.2f}", exp(b), -a)
~~~ 

*/