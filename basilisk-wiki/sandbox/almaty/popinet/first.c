/**
# Approximation of the first derivative
*/

#include "grid/cartesian1D.h"

int main() {
  init_grid(10);
  scalar f[];
  foreach()
    f[] = cos(2.*pi*x);
  boundary({f});
  
  scalar df[];
  foreach()
    df[] = (f[1,0] - f[0,0])/Delta;
  
  foreach()
    printf ("%g %g %g\n", x, f[], df[]);
  
  /**
  ~~~gnuplot f(x)
  set output 'f.png'
  df(x)=-2.*pi*sin(2.*pi*(x+1./40.))
  plot 'out' w lp, cos(2.*pi*x), 'out' using 1:3 w lp, df(x), 'out' u 1:($3-df($1)) w lp
  ~~~
  */
  
  free_grid();
}
