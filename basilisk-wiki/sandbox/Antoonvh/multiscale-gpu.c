/**
# Poisson problem on GPU

This poisson problem does not converge within the TOLERANCE after 100
MG cycles on a GPU grid, but is does on a CPU grid.
 */

double xc, yc, A, sigma;
double source_2D (double x, double y) {
  return A*4*(sq(x - xc) + sq(y - yc) - sq(sigma))*exp((-sq(x - xc) - sq(y - yc))/sq(sigma))/(sq(sq(sigma)));
}

double solution (double x, double y) {
  return A*exp((-sq(x - xc) - sq(y - yc))/sq(sigma));
}

//#include "grid/multigrid.h"
#define JACOBI 1
#include "poisson.h"

#define GN (3)

int main() {
  TOLERANCE = 1e-5;
  init_grid (256);
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
	b[] += source_2D(x, y);
	sol[] += solution (x,y);
      }
    }
  poisson (a, b);
  double err = 0, me = 0;
  foreach(reduction(+:err), reduction(max:me)) {
    double le = fabs (a[] - sol[]);
    err += sq(Delta)*le;
    if (le > me)
      me = le;
  }
  printf ("%d %g %g\n", N, err, me);
}