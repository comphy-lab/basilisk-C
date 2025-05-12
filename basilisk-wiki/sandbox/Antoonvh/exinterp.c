/**
# 6th-order accurate face-average to vertex interpolation

The explicit equation for the vertex values $f_0$

$$\frac{1}{4}f_{-1} + f_0 + \frac{1}{4}f_1 = \frac{3}{4} \left( F_{-1} + F_0 \right)$$

from a face-averaged values $F_{-1}$ and $F_0$ at the lhs and rhs of the vertex, respectively, yields a 4th order accurate interpolation.

Here we solve:

$$\frac{1}{3}f_{-1} + f_0 + \frac{1}{3}f_1 = \frac{29}{36} \left( F_{-1} + F_0 \right) + \frac{1}{36} \left( F_{-2} + F_1 \right) $$

for a 6th order accurate method;

~~~gnuplot
set xr [8: 2048]
set logscale xy 2
set ylabel 'L1 error'
set grid
set size square
plot 'out' u 1:3, 1e7*x**(-6)
~~~

 */

#include "grid/bitree.h"
#include "nsf2.h"

double func (double x, double y) {
  return exp(-sq(x));
}

scalar sc[];
vertex scalar sv[];
int main() {
  L0 = 20;
  X0 = -L0/2;
  for (int N = 16; N <= 1024; N *= 2) {
    init_grid (N);
    foreach()
      sc[] = Gauss6_y(x, 0, Delta, func); // average
    boundary ({sc});
    vertex scalar rhs[], temp[];
    int i = 0;
    foreach() {
      temp[] = 0;
       sv[] = ((-sc[-2] + 7*(sc[-1] + sc[]) - sc[1])/12.);
      rhs[] = 29./36.*(sc[] + sc[-1]) + 1./36.*(sc[1] + sc[-2]) ;
    }
    do {
      i++;
      foreach() 
	sv[] = rhs[] - 1./3.*(sv[-1] + sv[1]);
    } while (change (sv, temp) > 1e-12);
    double e = 0;
    foreach_vertex() 
      e += Delta*fabs(sv[] - func (x, 0));
    printf ("%d %d %g\n", N, i, e);
  }
  
  /**
## Vertex to face average

The inverse operator defined by,

$$ \frac{11}{38} F_{-1} + F_0 +\frac{11}{38} F_{1} = \frac{27}{38}\left( f_0 + f_1 \right) + \frac{3}{38} \left( f_{-1} + f_{2} \right) $$

is 6th-order accurate.

~~~gnuplot
plot 'log' u 1:3, 1e6*x**(-6)
~~~
   */

  for (int N = 16; N <= 1024; N *= 2) {
    init_grid (N);
    foreach_vertex()
      sv[] = func (x, 0);
    boundary ({sv});
    vertex scalar rhs[], temp[];
    int i = 0;
    foreach() {
      temp[] = 0;
      sc[] = (sv[] + sv[1])/2.;
      rhs[] = 27./38*(sv[] + sv[1]) + 3./38.*(sv[-1] + sv[2]);
    }
    do {
      i++;
      foreach() 
	sc[] = rhs[] - 11./38.*(sc[-1] + sc[1]);
    } while (change (sc, temp) > 1e-12);
    double e = 0;
    foreach() 
      e += Delta*fabs(sc[] - Gauss6_y (x, 0, Delta, func));
    fprintf (stderr, "%d %d %g\n", N, i, e);
  }
}
