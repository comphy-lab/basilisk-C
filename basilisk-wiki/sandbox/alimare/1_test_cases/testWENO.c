#define BGHOSTS 2

#include "utils.h"
#include "../weno5.h"

double sigma = 0.1;

int main() {
  double x0 = 0.5;

// We replace the constant value of n with a loop going through 16, 32, 64, 128 and 256.

  for (int n = 16; n <= 256; n *= 2) {
    init_grid (n);
    
    scalar f[], dfexact[];
    foreach() {
      f[] = erf((x-x0)/sigma);
      dfexact[] = 2./sigma/sqrt(M_PI)*exp(-pow((x-x0)/sigma,2));
    }
    boundary ({f, dfexact});
    
    vector dw5[];
    get_weno_derivative(f,dw5);
    foreach()
      dw5.x[]/=Delta;
    if (n == 32) {
      FILE * fp = fopen("solution.dat","w");
      foreach()
        fprintf (fp, "%g %g %g \n", x, dw5.x[],  dfexact[]);
      fclose(fp);
    }
      
    scalar e5[];
    foreach() {
      e5[] = dw5.x[] - dfexact[];
    }
    
    printf ("%d %g\n", n, normf(e5).max);
    
    free_grid();
  }
  exit(1);
}

/**

~~~gnuplot Numerical derivative approximation
set xlabel 'x'
set ylabel 'df/dx'
sigma=0.1
x0=0.5
f(x)=2./sigma/sqrt(pi)*exp(-((x-x0)/sigma)**2)
p "solution.dat" u 1:2 t 'Weno 5th' w p, \
  f(x) t 'Exact' w l
~~~ 

~~~gnuplot Convergence: The derivative converges with one order less than the approximation of the face values
set logscale
set xlabel '1/dx'
set ylabel 'normf(e).max'
set format y "10^%+03T"
fit [3:]a*x+b 'out' u (log($1)):(log($2)) via a,b
set xtics 8,2,512
set cbrange [1:2]
plot [8:512]'out' u 1:2 pt 7 t '5th-order WENO', \
            exp(b)*x**a t sprintf("%.0f/n^{%4.2f}", exp(b), -a)
~~~ 

*/ 