/**
#A test for the Poisson solver

This page concerns the convergence of the Poisson solver for $f(x,y)$,

$$\nabla^2 f = s.$$

Since [sine-rich
solutions](http://gerris.dalembert.upmc.fr/gerris/tests/tests/poisson.html) are
*not* suitable for testing, we design another test case where $s$ is
a dipole,

$$s = e^{-(x - 2)^2 - y^2} - e^{-(x + 2.)^2 - y^2} $$ 
 */
#include "grid/multigrid.h" 
#include "poisson.h"

scalar s[], f[];

#define SOURCE (exp(-sq(x - 2.) - sq(y)) - exp(-sq(x + 2.) - sq(y)))
/**
The cell-averaged source term is approximated with 64 points.
*/
double source (Point point) {
  double a = 0;
  foreach_child() foreach_child() foreach_child() a += SOURCE;
  return a/pow(pow(2., (double)dimension), 3.);
}
/**
# A reference solution

   Since it may be tricky to find an expression for the
   cell-averaged values of $f$, we compute a reference solution using a
   superior resolution grid. The solution on all relevant levels is
   then stored in an array.
 */

int maxlevel = 9;
int reference_level = 11;

double * sol;

void obtain_reference_solution () {
  long c = 0;
  int lr = 0;
  while (lr <= maxlevel)
    c += (1 << (dimension*lr++));
  sol = (double*)malloc(c*sizeof(double));
  init_grid (1 << reference_level);
  foreach() 
    s[] = source (point);
  poisson (f, s);
  restriction ({f});
  long ind = 0;
  for (int l = 0; l <= maxlevel; l++)
    foreach_level(l)
      sol[ind++] = f[];
}
/**
## The setup

Periodic conditions are set to exclude the error-inducing effects of
the approximate values of ghost cells. Furthermore we use a large
domain.
*/
int main() {
  periodic (bottom);
  periodic (left);
  TOLERANCE = 1e-7;
  L0 = 20.;
  X0 = Y0 = -L0/2;
  obtain_reference_solution ();
  /**
     We study the convergence starting from a 3-level grid down to the
`maxlevel`.
   */
  for (int l = 3; l <= maxlevel; l++) {
    init_grid (1 << l);
    foreach() 
      s[] = source (point);
    poisson (f, s);
    /**
       We skip to the relevant part of the array and compute a total
       error.
     */
    int lr = 0; long ind = 0;
    while (lr < l)
      ind += (1 << (dimension*lr++));
    double err = 0;
    foreach_level (l) 
      err += fabs(f[] - sol[ind++])*sq(Delta);
    printf ("%d\t%g\n", N, err);
  }
  free (sol);
}
  
/**
## Result

We obtain second order-accuracy

~~~gnuplot
set xr [4 : 1024]
set logscale xy 4
set size square
set grid 
set xlabel 'N'
set ylabel 'Total abs. error'
plot 'out' t 'data', 1000*x**-2 t 'second order'
~~~
 */ 
