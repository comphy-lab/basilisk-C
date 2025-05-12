/**
# The Helmholtz filter in 1D

We may filter a signal and add a boundary condition on the filtered
field. The filter also works on tree grids. 
 */

#include "grid/bitree.h"
#include "poisson.h"


int Helmholtz_filter (scalar phi, scalar phif, double Delta) {
  double sq_alpha = -sq(Delta/(2.*pi));
  const face vector alphaf[] = {sq_alpha, sq_alpha};
  return poisson (phif, phi, alphaf, unity).i;
}

scalar m[], mf[];
mf[left] = dirichlet (0);

int main() {
  init_grid (N);
  foreach()
    m[] = sq(sin(2*pi*x));
  boundary ({m});
  while (adapt_wavelet ({m}, (double[]){0.01}, 99).nf) {
    foreach()
      m[] = sq(sin(2*pi*x));
    boundary ({m});
  }
  Helmholtz_filter (m, mf, 0.5);
  foreach()
    printf ("%g\t%g\t%g\n", x, m[], mf[]);
}
/**
## The result

~~~gnuplot Great succes
plot 'out' u 1:2 w l t 'Original', \
     'out' u 1:3 w l t 'Filtered'
~~~
*/
