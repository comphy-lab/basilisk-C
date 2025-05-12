/** 
# Plotting a discrete derivative 

Given an error function compute the discrete derivative and 
compare it with the analytical value
*/

#include "grid/cartesian1D.h"
#include "utils.h"

int main() {
    
    int n = 100;
    double sigma = 0.1;
    double x0 = 0.3;
    init_grid (n);

    scalar f[], dfexact[];
    foreach() {
      f[] = erf((x-x0)/sigma);
      dfexact[] = 2./sigma/sqrt(M_PI)*exp(-pow((x-Delta/2.-x0)/sigma,2));
    }
  
    scalar df[];
    foreach()
      df[] = (f[] - f[-1])/Delta;

    foreach()
        printf ("%g %g %g \n", x-Delta/2., df[], dfexact[]);
}

/** After the program has finished we plot the function and the error

~~~gnuplot Numerical derivative approximation
set output 'derivative.png'
set xlabel 'x'
set ylabel 'df/dx'
p "out" u 1:2 t 'Numerical' w p pt 7, "out" u 1:3 t 'Exact' w l
~~~ 

~~~gnuplot Error of the numerical derivative approximation
set output 'error.png'
set xlabel 'x'
set ylabel 'err(df/dx)'
p "out" u 1:(abs($2-$3)) t 'Numerical' w lp pt 7
~~~ 

*/  