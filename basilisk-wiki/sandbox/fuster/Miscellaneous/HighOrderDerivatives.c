/**
# Discrete derivative 3rd-, 4rd and 5th-order 
*/

#include "grid/multigrid.h"
#include "utils.h"
#include "liblegendre.h"


/* Test function */
double sigma = 0.1;
double x0 = 0.5;

double fexact (double x) {

   return erf((x-x0)/sigma);

}


int main() {

  char fname[100];

  double xder[3] = {0.,0.,0.};

  for (int order = 3; order <= 5; order++) {

    sprintf(fname, "convergence.%i", order);
    FILE * fp0 = fopen(fname,"w");

    int r = 1 + (order > 4);
    LData MDer = Matrix_Derivate( xder, order, r );

    for (int n = 16; n <= 256; n *= 2) {
      init_grid (n);

      scalar fe[];
      vector dfe[];
      foreach() {
        fe[] = fexact(x);
        dfe.x[] = 2./sigma/sqrt(M_PI)*exp(-pow((x-x0)/sigma,2));
        dfe.y[] = 0.;
      }

      boundary((scalar *){fe, dfe});

      scalar favg[], dfeavgx[];
      Init_avg_variable ( fe, favg, order );
      Init_avg_variable ( dfe.x, dfeavgx, order );

      vector dw[];
      Variable_Gradient (favg, dw, &MDer);

      if (n == 32) {
        sprintf(fname, "solution.%i", order);
        FILE * fp = fopen(fname,"w");
        foreach()
          fprintf (fp, "%g %g \n", x, dw.x[]);
        fclose(fp);
      }

      scalar e[];
      foreach() 
        e[] = dw.x[] - dfeavgx[];

      fprintf (fp0, "%d %g\n", n, normf(e).max);

      free_grid();

    }

    deallocate_LData(&MDer);
    fclose(fp0);

  }

}

/**

~~~gnuplot Numerical derivative approximation
set xlabel 'x'
set ylabel 'df/dx'
sigma=0.1
x0=0.5
f(x)=2./sigma/sqrt(pi)*exp(-((x-x0)/sigma)**2)
p "solution.5" u 1:2 t 'deriv 5th' w p, \
  "solution.4" u 1:2 t 'deriv 4rd' w p, \
  "solution.3" u 1:2 t 'deriv 3rd' w p, \
  f(x) t 'Exact' w l
~~~ 

~~~gnuplot Convergence: The derivative converges with one order less than the approximation of the face values
set logscale
set xlabel '1/dx'
set ylabel 'normf(e).max'
fit [3:]a*x+b 'convergence.5' u (log($1)):(log($2)) via a,b
fit [3:]a1*x+b1 'convergence.4' u (log($1)):(log($2)) via a1,b1
fit [3:]a2*x+b2 'convergence.3' u (log($1)):(log($2)) via a2,b2
set xtics 8,2,512
set cbrange [1:2]
plot [8:512]'convergence.5' u 1:2 pt 7 t '5th-order deriv', \
            exp(b)*x**a t sprintf("%.0f/n^{%4.2f}", exp(b), -a), \
            'convergence.4' u 1:2 pt 7 t '4rd-order deriv', \
            exp(b1)*x**a1 t sprintf("%.0f/n^{%4.2f}", exp(b1), -a1), \
            'convergence.3' u 1:2 pt 7 t '3rd-order deriv', \
            exp(b2)*x**a2 t sprintf("%.0f/n^{%4.2f}", exp(b2), -a2)
~~~ 

*/