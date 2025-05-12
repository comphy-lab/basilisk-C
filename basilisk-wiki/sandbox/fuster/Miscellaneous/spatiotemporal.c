/** Spatio temporal solver */

#include "grid/multigrid.h"
#include "poisson.h"
#include "run.h"

scalar f[];
static double residual_xt (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  double maxres = 0.;
#if TREE
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*(a[] - a[-1])/Delta;
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*(a[] - a[0,-1])/Delta + (g.x[] - g.x[1])/Delta;
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*(a[] - a[0,-1])/Delta +
      ((alpha.x[1] + alpha.x[])*a[]
       - alpha.x[1]*a[1] - alpha.x[]*a[-1])/sq(Delta);
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif
  boundary (resl);
  return maxres;
}

static void relax_xt (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;

  foreach_level_or_leaf (l)
    a[] = (- sq(Delta)*b[] - lambda[]*Delta*a[0,-1] +
	   alpha.x[1]*a[1] + alpha.x[]*a[-1])
    /(- lambda[]*Delta + alpha.x[1] + alpha.x[]);
}

mgstats solver_xt (struct Poisson p)
{

  if (!p.alpha.x.i) {
    const vector alpha[] = {1.,1.,1.};
    p.alpha = alpha;
  }
  if (!p.lambda.i) {
    const scalar lambda[] = 1.;
    p.lambda = lambda;
  }

  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  face vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction ({alpha,lambda});

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ({a}, {b}, residual_xt, relax_xt, &p);

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}

f[left]    = dirichlet(1);
f[right]    = dirichlet(2);
@define dirichlettime(x) (x);
f[bottom]    = dirichlettime(0);



int main() {

    X0 = 0.;

    TOLERANCE = 1.e-6;

    for (int N =32; N <= 1024; N *= 2) {
      init_grid (N);
      run();
      free_grid();
    }
}

event output ( i = 0 ) {

  scalar r[], lambda[];
  face vector D[];
  foreach_face ()
    D.x[] = 1.;
  foreach() {
    lambda[] = -1.;
    r[] = -2.;
    f[] = 0.;
  }
  boundary({f});

  solver_xt (f, r, D, lambda);
  
  update_perf();

  double err = 0., errmax = 0.;
  foreach () {
    if (y > 1. - 1./N) {
      err   += fabs( f[] - ( 1. + 2*x - sq(x) ) )*Delta;
      errmax = max(errmax, fabs( f[] - ( 1. + 2*x - sq(x) ) ));
    }
  }

  double errtrans = 0.;
  foreach () {
    if (y + Delta/2. < 0.01 + 1.5/N && y + Delta/2. > 0.01 - 1.5/N) {
      if (x < 0.1) 
        errtrans = max(errtrans, fabs(erf(x/(2.*sqrt(y + Delta/2.))) -1. + f[]));
      if (y > 0.9) 
        errtrans = max(errtrans, fabs(erf((1.-x)/(2.*sqrt(y + Delta/2.))) - (2. - f[])/2.));
    }
  }

  fprintf(stderr, "%i %g %g %g %g %g \n", N, errtrans, errmax, err, perf.t, perf.speed);

  FILE * fp = fopen ("solution.ppm", "w");
  output_ppm(f, fp);
  fclose(fp);


}

/**
![Adaptive grid](spatiotemporal/solution.ppm)
*/

/**

~~~gnuplot error
unset y2label
unset y2tics
set xlabel 'N'
set ylabel 'errs'
set log xy
p "log" u 1:3 t 'steady' w lp, "log" u 1:2 t 'transient' w lp, '../spatiotemporal1/log' u 1:4 t 'Classical steady' w lp, '../spatiotemporal1/error.dat' u 1:2 t 'Classical transient' w p
~~~ 


*/
