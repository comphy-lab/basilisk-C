/* This is similar to gerris/test/poisson/circle */

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

void solve (int depth)
{
  size (1. [0]); // dimensionless space
  origin (-0.5, -0.5);
  int nrelax = 4;
  init_grid(1);

  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));

  foreach() {
    b[] =  -18.*sq(pi)*sin(3.*pi*x)*sin(3.*pi*y)*(x + y + 2.) +
      3.*pi*cos(3.*pi*x)*sin(3.*pi*y) +
      3.*pi*sin(3.*pi*x)*cos(3.*pi*y);
    a[] = 0.;
  }

  foreach_face()
    c.x[] = x + y + 2.;
  restriction ((scalar *){c});
  
  const scalar lambda[] = 0.;
  struct Poisson p = {0};
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
    foreach()
      if (fabs(res[]) > max)
	max = fabs(res[]);
    iter[i] = clock();
    maxres[i] = max;
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "residual %d %d %.1g\n", depth, i, maxres[i]);
    printf ("speed %d %d %g %g\n", depth, i, 
	    (iter[i] - start)/(double)CLOCKS_PER_SEC, maxres[i]);
  }

  double max = 0;
  foreach() {
    double e = a[] - solution(x, y);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, a[], b[], res[], e);
  }
  fprintf (stderr, "max error %d %g\n", depth, max);
}

int main (int argc, char ** argv)
{
  for (int depth = 7; depth <= 10; depth++)
    solve (depth);
}
