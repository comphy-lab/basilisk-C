/**
# Implicit (compact) ENO

For the development of a implicit (compact) WENO scheme, we first
investigate an implicit ENO scheme. Further, the "smoothness
indicator" is user knowledge.

We combine two compact 3-rd order  accurate one-sided (left and right) and a compact 4th order accurate centered schemes. 

~~~gnuplot Non smooth data
set xlabel 'x'
set ylabel 's'
plot 'out' u 1:2
~~~

~~~gnuplot Derivative
set xlabel 'x'
set ylabel 'ds/dx'
plot 'out' u 1:3
~~~

Well done `solve.h`!

*/
#define BGHOSTS 2
#include "grid/multigrid1D.h"
#include "solve.h"

#define CS (0.25*(ds[-1] + ds[1]) + ds[])
#define RCS (ds[] + 2*ds[1])
#define LCS (ds[] + 2*ds[-1])

#define RC (3./(4.*Delta)*(s[1 ] - s[-1 ]))
#define RRS ((-15./6.*s[0] + 2.*s[1] + s[2]/2.)/Delta)
#define LRS ((15./6.*s[0] - 2.*s[-1] - s[-2]/2.)/Delta)

int main() {
  L0 = 2;
  X0 = -1;
  TOLERANCE = 1e-9;
  periodic (left);
  init_grid (N);
  scalar s[], ds[];
  foreach() {
    ds[] = 0;
    s[] = x > 0 ? pi/2*x - pi/2 : (cos(pi*x/2));
  }

  // The "smoothness indicator" is apriori known
  solve(ds,
	((x > 0 && x < Delta) ? RCS : ((x < 0 && x > -Delta) ? LCS : CS)) ,
	((x > 0 && x < Delta) ? RRS : ((x < 0 && x > -Delta) ? LRS : RC)) );
  foreach()
    printf ("%g %g %g\n", x, s[], ds[]);
}
