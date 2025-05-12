/**
# A test for the Poisson solver on an adapted grid

Inspired by the centered Navier-Stokes solver, we aim to find the
cell-centered values of $f_x$, where $f$ statisfies a Poisson problem:

$$f_{xx} = s,$$

with $s$ a known source term. On this page we test with $s =
e^{-x^2}$, supplemented with boundary conditions for $f$.
*/
#include "grid/bitree.h"
#include "poisson.h"
scalar f[], s[];
f[left] = neumann (0.);
f[right] = dirichlet (0.);
#define SOURCE  (sqrt(pi)*(erf(x + Delta/2) - erf(x - Delta/2))/(2*Delta))
#define sol(x, c1, c2) ((c1*(x)) + c2 + sqrt(pi)*x*erf(x)/2. + exp(-sq(x))/2.)
#define DER ((1/Delta)*(sol((x + Delta/2.), c1, c2) - sol((x - Delta/2.),c1, c2)))
double c1, c2;
/**
We look at how the error convergens when we *naively* adapt the grid
based on $s$ with a decreasing limit on the so-called wavelet-based
discretization error estimate $\chi$ ($\zeta$).
 */
int main() {
  FILE * fp = fopen ("data", "w");
  L0 = 10.;
  X0 = -5.;
  /**
The constants $c1$ and $c2$ are determined by the boundary conditions
  */
  c1 = sqrt(pi)/2.;
  c2 = -sol (5, c1, 0);
  for (double zeta = 0.1; zeta > 0.0001; zeta /= 2){
    init_grid (1 << 3);
    do{
      foreach()
	s[] = SOURCE;
      boundary({s});
    }while (adapt_wavelet({s}, (double[]){zeta}, 19).nf);
    poisson (f, s);
    /**
We compute the cell-centered derivative (as is also done in
`centered.h`) and evaluate L1 error norm for the solution and the
derivative field.
     */
    double L1 = 0, L = 0;
    boundary ({f});
    foreach() {
      L1 += fabs((f[1] - f[-1])/(2*Delta) - DER)*Delta;
      L += fabs(f[] - sol(x ,c1, c2))*Delta;
    }
    fprintf (fp, "%ld\t%g\t%g\n", grid->n, L1, L);
  }      
  /**
We can study the convergence as a function of the used cells:

~~~gnuplot There appear to exist unfortunate grid structures
 set logscale xy
 set xlabel 'Cells'
 set ylabel 'Error'
 plot 'data' u 1:2 t 'Err. Derivative', 'data' u 1:3 t 'Err. Norm', 20*x**(-2) t 'second order convergence'
~~~

The results improve when $f$ is chosen to adapt upon rather than $s$. In the future, it will be shown how this bahaviour arises and how it may be mitigated. 
   */
} 