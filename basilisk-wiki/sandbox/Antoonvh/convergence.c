/**
# Convergence with adaptive refinement

The convergence of the $L_1$ and $L_\infty$ error norms for the
estimates of a face value and the cell centered first and second
derivatives in a scalar field with increasingly strict refinement
criteria is tested. We use second order accurate methods, which we
verify using the equidistant grid. The test function is a compact
Gaussian pulse, which is suitable for testing as it contains dominant
higher-order contributions.

First we define some helper functions,
 */
#include "grid/bitree.h"
#define FUNC (exp(-sq(x)))
#define FUNC1 (-2*x*FUNC)
#define FUNC2 ((4*sq(x) - 2)*FUNC)

#define D2SDX2 (-(s[-2] + s[2])/12. + 4.*(s[1] + s[-1])/3. - 5.*s[]/2.)

scalar s[];

/**
A function that evaluates the ability of a grid to reprent the
Gaussian and its derivatives
*/
int get_errors (double err[6]) {
  foreach()
    s[] = FUNC;
  boundary ({s});
  for (int i = 0; i < 6; i++)
    err[i] = 0;
  foreach_face() {
    double efl = fabs(FUNC - face_value(s, 0)); 
    if (err[3] < efl)
      err[3] = efl;
    err[0] += Delta*efl;
  }
  foreach() {
    double e1l = fabs(FUNC1 - (s[1] - s[-1])/(2*Delta));
    double e2l = fabs(FUNC2 - (s[1] - 2*s[] + s[-1])/sq(Delta));
    if (err[4] < e1l)
      err[4] = e1l;
    if (err[5] < e2l)
      err[5] = e2l;
    err[1] += Delta*e1l;
    err[2] += Delta*e2l;
  }
  return grid->tn;
}
/**
A function that prints the error data
 */
void print_errors (FILE * fp) {
  double err[6];
  fprintf (fp, "%d", get_errors (err));
  for (int i = 0; i < 6; i++)
    fprintf (fp, " %g", err[i]);
  fputc ('\n', fp);
}
/**
First, a simple convergence test.
*/
int main() {
  L0 = 20;
  X0 = -L0/2.;
  for (int l = 3; l <= 10; l++) {
    init_grid (1 << l);
    print_errors (stdout);
  }
  /**
The result is:
     
~~~gnuplot check 2nd order scheme test
set xr [4:2048]
set xlabel 'N'
set ylabel 'error'
set logscale y 
set logscale x 2
set grid
plot 'out' u 1:2, '' u 1:3, '' u 1:4, '' u 1:5, '' u 1:6, '' u 1:7,  200*x**(-2) lw 2 t 'second order'
~~~

Now the grid is adaptively refined with an increasingly stringent criterion. 
   */
  for (double zeta = .1; zeta > 5e-5; zeta /= 2) {
    init_grid (8);
    do {
      foreach()
	s[] = FUNC;
      boundary({s});
    } while (adapt_wavelet({s}, (double[]){zeta}, 99).nf);
    print_errors (stderr);
  }
  /**
We check the convergence for each of the 3 tests

~~~gnuplot Face interpolator error
set xr [16:512]
set xlabel 'cells'
set ylabel 'error'
set logscale y 
set logscale x 2
set grid
plot 'log' u 1:2 t 'L_1 adaptive', '' u 1:5 t 'L_{inf} adapte', 200*x**(-2) lw 2 t 'second order'
~~~

Good...

~~~gnuplot first derivative error
set xr [16:512]
set xlabel 'cells'
set ylabel 'error'
set logscale y 
set logscale x 2
set grid
plot 'log' u 1:3 t 'L_1 adaptive', '' u 1:6 t 'L_{inf} adapte',\
 200*x**(-2) lw 2 t 'second order'
~~~

Ok....

~~~gnuplot second derivative error
set xr [16:512]
set xlabel 'cells'
set ylabel 'error'
set logscale y 
set logscale x 2
set grid
set key bottom left
plot 'log' u 1:4 t 'L_1 adaptive', '' u 1:7 t 'L_{inf} adapte',\
 200*x**(-2) lw 2 t 'second order',  10*x**(-1) lw 1 t 'First order'
~~~

## Wut?

Consider the stencil,

$$\frac{\mathrm{d}^2s}{\mathrm{d}x^2} \approx \frac{s[-1] - 2s[0]
+ s[1]}{\Delta^2} + c_2 \frac{\mathrm{d}^4s}{\mathrm{d}x^4}   \Delta^{-2}$$

is second order accurate. But the neighbor values need to be estimated
at resolution boundaries where they pollute the approximation:

$$\left. \frac{\mathrm{d}^2s}{\mathrm{d}x^2}
\right\vert_{\mathrm{res.\ bound.}} \approx \frac{s[-1] - 2s[0] + s[1]
+ c\frac{\mathrm{d}^2s}{\mathrm{d}x^2}\Delta^2}{\Delta^2} \approx
\frac{s[-1] - 2s[0] + s[1]}{\Delta^2} + \frac{\zeta}{\Delta^2},$$

Where $c$ is some contant and $\zeta$ is the refinement
criterion. This extra term scales with $\Delta^{-2}$, meaning that is
gets bigger with decreasing $\Delta$! However, the local grid-cell
($\Delta$) size is a function of $\zeta$. The `adapt_wavelet`
algorithm tries to maintain a certain wavelet-based error estimation
($\xi$):

$$\xi \approx c\frac{\mathrm{d}^2s}{\mathrm{d}x^2}\Delta^2 \approx
\zeta = \mathrm{constant},$$

and thus,

$$\Delta(\zeta) =
\sqrt{\frac{\zeta}{c\frac{\mathrm{d}^2s}{\mathrm{d}x^2}}}.$$

Meaning for the estimate near resolution boundaries:

$$\left. \frac{\mathrm{d}^2s}{\mathrm{d}x^2}\right\vert_{\mathrm{res. bound.}} \approx \frac{s[-1] - 2s[0]
+ s[1]}{\Delta^2} + c\frac{\mathrm{d}^2s}{\mathrm{d}x^2},$$

which is zero-th order accurate.  

We can have a look at the error distribution:
  */
  boundary ({s});
  FILE * fp = fopen ("level", "w");
  foreach()
    fprintf (fp, "%g %d %g\n", x, level,
	     fabs(FUNC2 - (s[1] - 2*s[] + s[-1])/sq(Delta)));
  fclose (fp);
  /**

~~~gnuplot The errors appear to coincide with resolution boundaries.
reset
set xr [-4:4]
set samples 400
set grid
plot (4*x**2 - 2)*exp(-x**2) t '2nd derivative', \
 4*(4*x**4 - 12*x**2 + 3)*exp(-x**2) t '4th' ,\
 'level', 'level' u 1:(25*$3) t 'Error'
~~~

So if we would refine all cells by one level, we will see proper
convergence. The problem is that that grid structure is not found by
`adapt_wavelet()`.

The analysis is repeated with a 4th-order accurate approximation
  */
  boundary ({s});
  FILE * fp2 = fopen ("4th", "w");
  foreach()
    fprintf (fp2, "%g %d %g\n", x, level,
	     fabs(FUNC2 - D2SDX2/sq(Delta)));
  fclose (fp2);
  /**
~~~gnuplot
reset
set logscale y
set xr [-4:4]
set grid
set ylabel 'error' 
plot 'level' u 1:3, '4th' u 1:3 t 'Error'
~~~

The issue is with the treatment of the resolution boundary, not the
approximation itself.

## Is this an issue?

The idea that the refinement criterion can be tuned to obtain more
accurate approximations is not true. If you never believed that, there
is no issue.

## is there a solution?

Yes! If we make a proper analysis of the error, it is possible to
achieve hyper-convergence for certain error-norms. That maybe
worthwhile.
*/
}
