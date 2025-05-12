/**
# 4th-order compact central/upwind finite difference

Solving,

$$\frac{17}{14} f'_{-1} + f'_{0} + \frac{-1}{14}f' = \frac{\frac{27}{14}f_{0} + \frac{-12}{7}f_{-1} + \frac{-3}{14}f_{-2}}{\Delta},$$

yields a 4th-order accurate scheme for $f'$. So does,

$$f'_{-2} + 4 f'_{-1} + f'_0 = \frac{3f_{0} -3f_{-2}}{\Delta}.$$


~~~gnuplot
set logscale xy
set xr [8 : 1024]
set grid
set size square
set xlabel 'N'
set ylabel 'L1 error'
plot 'out't 'data', 1e5*x**(-4) t '4th order'
~~~
*/
#include "grid/bitree.h"

void upwcd (scalar s, scalar da) {
  double alpha = 1., beta = 4; // Larger than 1!
  double a = -3., b = 0., c = 3;
  scalar rhs[];
  foreach() {
    rhs[] = (c*s[] + b*s[-1] + a*s[-2])/Delta;
    da[] = (s[1] - s[-1])/(2*Delta);
  }
  for (int i = 0; i < depth(); i++) {
    foreach() {
      double n = rhs[] - alpha*da[-2] - beta*da[-1];
      da[] = (5*da[] + n)/6.;
    }
    
    boundary ({da});
  }
}
int compact_iters = 10;
void compact_first_derivative (scalar * sl, vector * dsl) {
  vector * rhsl;
  rhsl = (vector*)list_clone ((scalar*)dsl);
  double alphaP = 1./4., aP = 1.5;
  foreach() {
    scalar s;
    vector ds, rhs;
    for  (s, ds, rhs in sl, dsl, rhsl) {
      foreach_dimension() {
	ds.x[] = 0;//((s[1] - s[-1])/(2*Delta));
	rhs.x[] = aP*(s[1] - s[-1])/(2.*Delta);
      }
    }
  }
  boundary((scalar*)dsl);
  for (int it = 0; it < compact_iters; it++) {
    foreach() {
      vector ds, rhs;
      for  (ds, rhs in dsl, rhsl) {
	foreach_dimension()	
	  ds.x[] = rhs.x[] - alphaP*(ds.x[-1] + ds.x[1]);
      }
    }
    boundary ((scalar*)dsl);
  }
  delete((scalar*)rhsl);
  free (rhsl); 
  rhsl = NULL;
}


#define FUNC (exp(-sq(x)))
#define SOL (-2*x*exp(-sq(x)))

int main() {
  periodic (left);
  L0 = 25;
  X0 = -L0/2;
  for (int N = 16; N <= 512; N *= 2) {
    init_grid (N);
    scalar s[], ds[];
    foreach()
      s[] = FUNC;
    boundary ({s});
    compact_first_derivative ({s}, (vector*){ds});
    double e = 0;
    foreach()
      e += Delta*fabs(ds[] - SOL);
    printf ("%d %g\n", N, e);
  }
}
