/**
# Convergence for the first derivative of a function

We start with the code from the [previous exercise](first.c) and modify it
to compute the *convergence rate* of the approximation. */

#include "grid/cartesian1D.h"
#include "utils.h"

int main() {
  /** We replace the constant value of *n* with a loop going through
  10, 20, 40, 80 and 160. */
  for (int n = 10; n <= 160; n *= 2) {
    /** The *body* of the loop is just copied from the first exercise... */
    init_grid (n);

    scalar f[];
    double k = 2.*pi;
    foreach()
      f[] = cos(k*x);

    scalar df[];
    foreach()
      df[] = (f[] - f[-1])/Delta;

    /** ... except that we don't want a graph of df(x) for each value of *n*. What we
    would like is a graph of a *measure* of the error as a function of *n*. We
    can do that by computing a *norm* of the error, for example
    $$
    \max_x(|df(x) + k \sin(k x)|)
    $$
    the maximum error between the numerical and analytical value of the derivative. 
    
    We first compute a new field *e* to hold the error for each grid point. */

    scalar e[];
    foreach()
      e[] = df[] + k*sin(k*(x - Delta/2.));

    /** We can then call the function *normf()* (defined in the header file *utils.h*
    which we included above), to get the maximum error (as well as other norms).
    
    We use the *printf()* function to output *n* and the error (the '%d' formatting 
    character means that *n* is an integer rather than a floating-point number). */
    
    printf ("%d %g\n", n, normf(e).max);
  }
}

/** After the program has finished we plot the error as a function of *n*

~~~gnuplot
set logscale
set xtics 5,2,320
set xlabel 'n'
set ylabel 'normf(e).max'
fit a*x+b 'out' u (log($1)):(log($2)) via a,b
plot [5:320]'out' pt 7 t '', exp(b)*x**a t sprintf("%.0f/n^{%4.2f}", exp(b), -a)
~~~ 

We use a log-log scale, which clearly shows that the error follows a power law
(i.e. a straight line on a log-log graph), with an exponent close to -2 as revealed 
by the gnuplot fit (have a look at the [raw page source](/_showraw/sandbox/course/first2.c)
for the gnuplot commands). This means that our estimate of the first derivative converges at
a *second-order* rate with spatial resolution. In other words, doubling the resolution will
divide the error by four. */