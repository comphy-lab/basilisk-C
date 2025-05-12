/**
# Designing a multiscale problem for the Poisson solver.

For the design of a relevant multiscale-test case that can help to
assess tree-refinement indicators. Herein we study the multigrid
convergence of academic Poisson problem with source:

![The Poisson source on grids with increasing resolutions](multiscale/b.mp4)

~~~gnuplot Convergence is as expected.
set size square
set xr [200:5000]
set yr [0.001:1]
set logscale x 2
set logscale y
set xlabel 'N'
set ylabel 'Error'
set grid
plot 'out' u 1:2 t 'L1 error', '' u 1:3 t 'Max error', 2e4*x**(-2) t 'second order
~~~

 */

// Source
double xc, yc, A, sigma;
double source_2D (double x, double y) {
  return A*4*(sq(x - xc) + sq(y - yc) - sq(sigma))*exp((-sq(x - xc) - sq(y - yc))/sq(sigma))/(sq(sq(sigma)));
}

double solution (double x, double y) {
  return A*exp((-sq(x - xc) - sq(y - yc))/sq(sigma));
}

#include "grid/multigrid.h"
#include "poisson.h"
#include "utils.h"
#include "higher-order.h"

#define GN (3)

int main() {
  TOLERANCE = 1e-5;
  for (N = 256; N <= 4096; N *= 2) {
    init_grid (N);
    scalar a[], b[], sol[];
    foreach_dimension() {
      a[left] = dirichlet (0);
      a[right] = dirichlet (0);
    }
    L0 = 40;
    double sigmas[GN] = {.2, .5, 1};
    double As[GN] = {1,2,4};
    
    foreach() 
      a[] = b[] = sol[] = 0;
    
    for (int i = 0; i < GN; i++)
      for (int j = 0; j < GN; j++) {
	xc = X0 + (i + 1)*L0/(GN + 1);
	yc = Y0 + (j + 1)*L0/(GN + 1);
	sigma = sigmas[i];
	A = As[j];
	foreach() {
	  b[] += Gauss6 (x, y, Delta, source_2D);
	  sol[] += Gauss6 (x, y, Delta, solution);
	}
      }
    output_ppm (b, file = "b.mp4", n = 400, opt = "-r 2");
    poisson (a, b);
    output_ppm (a, file = "a.mp4", n = 400,
		opt = "-r 2", min = -2, max = 2);
    output_ppm (sol, file = "sol.mp4", n = 400,
		opt = "-r 2", min = -2, max = 2);
    double err = 0, me = 0;
    foreach(reduction(+:err) reduction(max:me)) {
      double le = fabs (a[] - sol[]);
      err += sq(Delta)*le;
      if (le > me)
	me = le;
    }
    printf ("%d %g %g\n", N, err, me);
  }
}
