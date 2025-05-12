/**
# 2 + 2 = 4

We combine two second-order accurate solutions to a Poisson equation
to get a 4th-order accurate one for the gradients.

![The source term to the test Poisson problem](corrected_poisson/source.png)

![The solution to the test problem](corrected_poisson/solution.png)

~~~gnuplot Error in the gradients
set logscale x 2
set logscale y
set xr [16:2048]
set yr [1e-5:10]
set size square
set grid
plot 'out' u 1:2 t 'face gradient', '' u 1:3 t 'corrected gradient', 1e7*x**(-4) t 'fourth order', 2e3*x**(-2) t 'second order'
~~~

 */
#define BGHOSTS 2
#include "grid/multigrid.h"
#include "poisson.h"
#include "higher-order.h"
#include "utils.h"

double solution (double x, double y) {
  return exp(-sq(2*sq(x) + sq(y) - 2));
}

double source (double x, double y) {
  return 4*solution(x, y)*(6 + 64*pow(x, 6) + 11*sq(y) -
			   16*pow(y, 4) + 4*pow(y, 6) +
			   16*pow(x, 4)*(5*sq(y) - 8) +
			   sq(x)*(50 - 96*sq(y) + 32*pow(y, 4)));
}

double grad_x (double x, double y) {
  return solution(x, y)*-8*x*(2*sq(x) + sq(y) - 2);
}

double grad_y (double x, double y) {
  return solution(x, y)*-4*y*(2*sq(x) + sq(y) - 2);
}

double error (face vector g) {
  double e = 0;
  foreach_face(reduction(+:e))
    e += sq(Delta)*fabs(g.x[] - Gauss6_x (x, y, Delta, grad_x));
  return e;
}

scalar a[], b[];

int main() {
  L0 = 12;
  X0 = Y0 = -L0/2.;
  TOLERANCE = 1e-9;
  for (N = 32; N <= 1024; N *= 2) {
    init_grid (N);
    a[left] = dirichlet (solution(x, y));
    foreach() {
      a[] = 0;
      b[] = Gauss6(x, y, Delta, source);
    }
    poisson (a, b);
    if (N == 512) {
      output_ppm (b, file = "source.png", n = 400);
      output_ppm (a, file = "solution.png", min = -1, max = 1,
		  map = blue_white_red, n = 400);
    }
    // Compute gradient field and its error
    face vector ga[];
    foreach_face()
      ga.x[] = face_gradient_x(a, 0);
    double err1 = error(ga);
    // Compute and apply gradient correction:
    vector gac[];
    foreach() {
      foreach_dimension()
	gac.x[] = (a[1] - a[-1])/(2*Delta);
    }
    // 4th-order pseudo circulation
    foreach()
      b[] = -1/6.*(gac.y[2] - gac.y[-2] - gac.x[0,2] + gac.x[0,-2])/(4*Delta);
    poisson (a, b);
    coord f = {-1.,1.};
    foreach_face()
      ga.x[] -= f.x*(a[0,1] + a[-1,1] - a[0,-1] - a[-1,-1])/(4.*Delta);
    printf ("%d %g %g\n", N, err1, error(ga));
  }
  
}
