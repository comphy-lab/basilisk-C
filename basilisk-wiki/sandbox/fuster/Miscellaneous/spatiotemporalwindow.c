/** Spatio temporal solver */

#include "poisson.h"
#include "run.h"

scalar f[];

@define dirichlettime(x) (x);
f[bottom]    = dirichlettime(0);
f[left]    = dirichlet(1);
f[right]   = dirichlet(2);

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

scalar r[], lambda[];
face vector D[];


int main() {

    int N = 4; 
    init_grid (N);

    run();

}

event init (i = 0 ) { 

    foreach_face ()
      D.x[] = 1.; 
    foreach() {
      lambda[] = -1.;
      r[] = -2.;
      f[] = 0.; 
    }   

}

event takestep (i++; t < 1.) {

  dt = dtnext (0.1);

    scalar frefine[];

      foreach () {
          frefine[] = f[]*exp(-sq(y-t)/sq(0.5*dt));
      }   
      adapt_wavelet ({frefine}, (double []){1e-2}, maxlevel = 8, minlevel = 2); 
      solver_xt (f, r, D, lambda);

}


event output (i++; t < 1.) {

      double err = 0., errmax = 0.;
      foreach () {
        if (y > 0.95) {
          err   += fabs( f[] - ( 1. + 2*x - sq(x) ) );
          errmax = max(errmax, fabs( f[] - ( 1. + 2*x - sq(x) ) ));
        }
      }


      fprintf(stderr, "%g %ld %g %g %g \n", t, grid->tn, perf.t, err, errmax);

      scalar l[];
      foreach () 
        l[] = level;

      output_ppm (l, file = "grid.mp4", n = 512);

        foreach () {
          if (y > 0.95)
            printf("%g %g %g \n", x, t, f[]);
        }

}
/**
 ![Evolution of the level of refinement](spatiotemporalwindow/grid.mp4)
 */

/**

~~~gnuplot Exact steady state solution
set xlabel 'x'
set cblabel 't'
set ylabel 'f'
set key left
p "out" u 1:3:2 not w p palette pt 7, 1 - x**2 + 2.*x t 'EXACT' w l lw 3 lc 0
~~~ 

~~~gnuplot number of cells
set xlabel 't'
set cblabel 'cells'
p "log" u 1:2 not w l
~~~ 

~~~gnuplot error vs. tcpu
set log xy
set xlabel 'tcpu'
set ylabel 'errs'
set key below
p "log" u 3:5 t 'current' w lp, "../spatiotemporal/log" u 5:4 t 'spatiotemporal' w lp
~~~ 

*/ 
