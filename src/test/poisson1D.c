#include "grid/multigrid1D.h"
#include "utils.h"
#include "poisson.h"

scalar a[], b[];

double solution (double x)
{
  return sin(3.*pi*x);
}

/* Dirichlet condition on all boundaries */
a[right]  = dirichlet (solution(x));
a[left]   = dirichlet (solution(x));

int main (int argc, char ** argv)
{
  size (1.[0]);
  origin (-0.5, -0.5);

  for (int depth = 5; depth <= 8; depth++) {
    init_grid (1 << depth);

    foreach() {
      b[] = - sq(3.*pi)*sin (3.*pi*x);
      a[] = 0.;
    }

    TOLERANCE = 1e-4;
    poisson (a, b);

    double max = 0;
    foreach() {
      double e = a[] - solution(x);
      if (fabs(e) > max) max = fabs(e);
      printf ("%g %g %g %g\n", x, a[], b[], e);
    }
    fprintf (stderr, "%d %g\n", depth, max);
  }
}

/** After the program has finished we plot the error as a function of *n*

~~~gnuplot
set logscale y
set xlabel 'level'
set ylabel 'normf(e).max'
fit a*x+b 'log' u 1:(log($2)) via a,b
plot [4:9]'log' pt 7 t '', exp(a*x+b) t sprintf("%.0f/n^{%4.2f}", exp(b), -a/log(2.))
~~~ 
*/
