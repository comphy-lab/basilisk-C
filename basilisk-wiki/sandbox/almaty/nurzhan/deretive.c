/**
# Approximation of the first derivative
*/

#include "grid/cartesian1D.h"

int main() {
  init_grid(10);
  scalar f[];
  foreach()
    f[] = sin(2.*pi*x);
  foreach()
    printf ("%g %g\n", x, f[]);
  
  /**
  ~~~gnuplot f(x)
  set output 'f.png'
  plot 'out'
  ~~~
  */
  
  free_grid();
}