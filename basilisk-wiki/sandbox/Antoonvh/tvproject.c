/**
## Projection of a vertex vector field

This is part of the "Hogde Decomposition" [test
case](/sandbox/jmf/Hodge/hodge.c) by Jose Maria Fullana.

~~~gnuplot Original
set size square
plot 'out' u 1:2:($3/10):($4/10) with vect t "original"
~~~

~~~gnuplot Projected
plot 'log' u 1:2:($3/10):($4/10) with vect t "projected"
~~~
 */

#include "nodal-poisson.h"
#include "utils.h"

vector v[];
scalar p[];

int main() {
  L0 = 4;
  X0 = Y0 = -L0/2;
  foreach_dimension()
    periodic(left);
  init_grid (32);
  p.prolongation = refine_vert;
  p.restriction = restriction_vert;
  foreach_dimension() {
    v.x.prolongation = refine_vert;
    v.x.restriction = restriction_vert;
  }
  foreach_vert() {
    p[] = 0;
    v.x[] = sin(pi*x)*cos(pi*y) + cos(pi*x)*sin(pi*y);
    v.y[] = cos(pi*x)*sin(pi*y) - sin(pi*x)*cos(pi*y);
  }
  refine ((fabs (v.x[]) > 0.5 || fabs(v.y[]) > 0.5) && level < 6);
  foreach_vert() 
    printf ("%g %g %g %g\n", x, y, v.x[], v.y[]);

  vproject (v, p, 1e-3);
  
  foreach_vert() 
    fprintf (stderr, "%g %g %g %g\n", x, y, v.x[], v.y[]);
}
