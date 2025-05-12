/**
# 8th-order accurate finite-differencing using a 5-point stencil

Solving the system 

$$\beta f'_{-2} + \alpha f'_{-1} + f'+\alpha f'_{1} + \beta f'_{2} =
a\frac{f_{1} - f_{-1}}{2 \Delta }+ b\frac{f_{2} - f_{-2}}{4\Delta},$$

with coefficients,

$$\alpha = \frac{4}{9},$$
$$\beta = \frac{1}{36},$$
$$a = \frac{40}{27},$$
$$b = \frac{25}{54},$$

yields an 8-th order accurate estimation for $f'$ (Lele, 1991).

~~~gnuplot Convergence is OK
set logscale xy 
set grid 
set size square 
set xr [4:1024] 
set xlabel 'N'
set ylabel 'error' 
plot 'out' u 1:2 t 'L1',\
     '' u 1:3 t 'Max' ,\
 1e10*x**-8 t '8th order' 
~~~

 ~~~gnuplot Convergence is OK
reset
set logscale x 
set grid 
set size square 
set xr [4:1024] 
set yr [0:9]
set xlabel 'N'
set ylabel 'Iterations' 
set key top left
plot 'out' u 1:4 t 'Cycles', '' u 1:5 t 'sweeps'
~~~
*/
#define BGHOSTS 2
#include "grid/multigrid1D.h"
#include "poisson.h"

#define DF1(f) ((f[1] - f[-1])/(2*Delta))
#define DF2(f) ((f[2] - f[-2])/(4*Delta))
double alphaP = 4./9., betaP = 1./36.;
double aP = 40./27., bP = 25./54.;

double residual_pade (scalar * al, scalar * bl,
		      scalar * resl, void * data) {
  scalar a = al[0], b = bl[0], res = resl[0];
  double mr = 0;
  foreach(reduction (+:mr)) {
    res[] = (a[] + alphaP*(a[1] + a[-1]) + betaP*(a[2] + a[-2]) -
	     (aP*DF1(b) + bP*DF2(b)));
    if (fabs(res[]) > mr)
      mr = fabs(res[]);
  }
  return mr;
}

void relax_pade (scalar * al, scalar * resl,
		 int depth, void * data) {
  scalar a = al[0], res = resl[0];
  foreach_level(depth) {
    a[] = -(res[] + alphaP*(a[-1] + a[1]) + betaP*(a[-2] + a[2]));
  }
}

/**
## Set-up a test case 

Differentiating a Gaussian pulse.
 */

#define F(x) (exp(-sq(x)))
#define DF(x) (-2.*x*exp(-sq(x)))

scalar f[], dfdx[];
int main() {
  L0 = 20;
  X0 = -L0/2.;
  for (N = 8; N <= 512; N *= 2) {
    init_grid(N);
    foreach() 
      f[] = F(x);
    boundary ({f});
    foreach()
      dfdx[] = 4./3.*DF1(f) - DF2(f)/3.; //Guess
    boundary ({dfdx}); // BC for the derivative...
    mgstats mgs = mg_solve ({dfdx}, {f},
			    residual_pade, relax_pade);
    double e = 0, em = 0;
    foreach() {
      double el = fabs(DF(x) - dfdx[]);  
      e += Delta*el;
      if (el > em)
	em = el;
    }
    printf ("%d %g %g %d %d\n", N, e, em, mgs.i, mgs.nrelax);
    TOLERANCE = 1e-3*em;
  }
}



/**
## Note

A 6th order accurate approximation for $f'$ can be obtained with the
coefficients.

$$a = \frac{2}{3}(\alpha + 2),$$
$$b = \frac{1}{3}(4\alpha - 1),$$
$$\alpha = \frac{1}{3},$$
$$\beta = 0.$$

Which has a 3-point relaxation stencil.

Futher, a 4th order method is achievable on a 3-point stencil.

## Reference

S.K. Lele, *Compact Finite Difference Schemes with Spectral-like Resolution
* (1991), in JCP.
*/
