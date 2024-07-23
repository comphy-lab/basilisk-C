// Poisson equation with variable coefficients

#include "utils.h"
#include "poisson.h"

scalar a[], b[], res[], dp[];
face vector c[];

double solution (double x, double y)
{
  return sin(3.*pi*x)*sin(3.*pi*y);
}

/* Dirichlet condition on all boundaries */
a[right]  = dirichlet (solution(x, y));
a[left]   = dirichlet (solution(x, y));
a[top]    = dirichlet (solution(x, y));
a[bottom] = dirichlet (solution(x, y));

/* homogeneous conditions for dp */
dp[right]  = dirichlet(0);
dp[left]   = dirichlet(0);
dp[top]    = dirichlet(0);
dp[bottom] = dirichlet(0);

int main (int argc, char ** argv)
{
  size (1. [0]);
  origin (-0.5, -0.5);
  int depth = argc < 2 ? 9 : atoi(argv[1]), nrelax = 4;
  init_grid(1 << depth);

  foreach() {
    b[] = -18.*sq(pi)*sin(3.*pi*x)*sin(3.*pi*y)*(x + y + 2.) +
      3.*pi*cos(3.*pi*x)*sin(3.*pi*y) +
      3.*pi*sin(3.*pi*x)*cos(3.*pi*y);
    a[] = 0.;
  }

  foreach_face()
    c.x[] = x + y + 2.;
  restriction ((scalar *){c});
  const scalar lambda[] = 0.;
  struct Poisson p;
  p.a = a; p.b = b; p.alpha = c; p.lambda = lambda;
  #define NITER 13
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  scalar * lres = {res};
  residual ({a}, {b}, lres, &p);
  for (int i = 0; i < NITER; i++) {
    mg_cycle ({a}, lres, {dp}, relax, &p, nrelax, 0, depth());
    residual ({a}, {b}, lres, &p);
    double max = 0.;
    foreach(reduction(max:max))
      if (fabs(res[]) > max)
	max = fabs(res[]);
    iter[i] = clock();
    maxres[i] = max;
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "%d %.2g\n", i, maxres[i]);
    printf ("%d %g %g\n", i, (iter[i] - start)/(double)CLOCKS_PER_SEC, 
    	    maxres[i]);
  }
  double max = 0;
  foreach(reduction(max:max)) {
    double e = a[] - solution(x, y);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, a[], b[], res[], e);
  }
  fprintf (stderr, "# max error %g\n", max);
}
