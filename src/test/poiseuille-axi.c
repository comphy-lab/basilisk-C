/**
# Axisymmetric Poiseuille flow

Same as [poiseuille-periodic.c]() but axisymmetric. */

#include "axi.h"
#include "navier-stokes/centered.h"

int main() {
  periodic (right);

  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */
  
  size (1. [0]);
  DT = HUGE [0];

  TOLERANCE = 1e-6;

  u.t[top] = dirichlet(0);

  for (N = 8; N <= 64; N *= 2)
    run(); 
}

scalar un[];

event init (t = 0) {
  a[] = {1.,0.};
  mu = fm;
  foreach()
    un[] = u.x[];
}

event logfile (t += 0.1; i <= 100) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-6)
    return 0; /* stop */
}

event profile (t = end) {
  printf ("\n");
  foreach()
    printf ("%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
  scalar e[];
  foreach()
    e[] = u.x[] - 0.25*(1. - y*y);
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}

/**
~~~gnuplot Second-order convergence is obtained as expected.
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2
set xlabel 'Resolution'
set ylabel 'Error'
set logscale
set xrange [4:128]
set cbrange [1:2]
set xtics 4,2,128
set grid ytics
set yrange [1e-5:]
plot 'log' u 1:4 t 'max', exp(f(log(x))) t ftitle(a,b), \
     'log' u 1:2 t 'norm1', exp(f2(log(x))) t ftitle(a2,b2)
~~~
*/
