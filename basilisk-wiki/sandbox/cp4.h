/**
# A compact 4th-order Finite volume Poisson solver

Using only a 3-point stencil.

We need a few higher-order methods
*/
#include "higher-order.h"

/**
The mg-cycle prolongation is conservative and 3rd-order accurate
 */

static inline double mul_3pt (Point point, scalar s) {
  return mul_xpt (point, s, 3);
}
#define bilinear mul_3pt
#include "solve.h"
#undef bilinear

/**
On trees, the face-gradient needs to be determined using an implicit scheme.

$$\frac{1}{10}f'_{i-1\frac{1}{2}} + f'_{i-\frac{1}{2}} +
\frac{1}{10}f'_{i+\frac{1}{2}} = \frac{6}{5\Delta}\left( f_i - f_{i-1}
\right)$$
 */
static double compact_der (face vector da, scalar a, double TOLER) {
  double maxres = HUGE;
  while (maxres > TOLER) {
    maxres = -1;
    foreach_face(reduction (max:maxres)) {
      double resa = 6.*(a[] - a[-1])/(5.*Delta) -
	(da.x[-1] + da.x[1])/10. - da.x[];
      if (fabs(resa) > maxres)
	maxres = fabs(resa);
      da.x[] += resa; 
    }
    boundary ((scalar*){da});
  }
  return maxres;
}
/**
The relaxation and residual function
*/
static void relax_compact (scalar * da, scalar * res,
			    int l, void * data) {
  scalar a = da[0], b = res[0];
#if JACOBI
  scalar c[];
#else
  scalar c = a;
#endif
  foreach_level_or_leaf (l) {
    double n = -5*sq(Delta)*b[]/6., d = 2*dimension;
    n += 4./5.*(a[1] +  a[-1] + a[0,1] + a[0,-1]) +
      1./5.*(a[1,1] + a[1,-1] + a[-1,-1] + a[-1,1]);
    c[] = n/d;
  }
  
  /**
     For weighted Jacobi we under-relax with a weight of 2/3. */
  
#if JACOBI
  foreach_level_or_leaf (l)
    a[] = (a[] + 2.*c[])/3.;
#endif
}

double _tol_prev = -1;
static double residual_compact (scalar * al, scalar * bl,
				scalar * resl, void * data) {
  scalar a = al[0], b = bl[0], res = resl[0];
  int i = 1;
  double maxres = -1;
#if TREE
  face vector da;
  foreach_dimension()
    da.x = bl[i++];
  boundary ((scalar*){da});
  /**
The tolerance for the face gradients scales with the residual for the
Poisson problem of the previous iteration.
   */
  if (_tol_prev < 0)
    _tol_prev = TOLERANCE;
  compact_der (da, a, _tol_prev/10.);
  foreach (reduction(max:maxres)) {
    res[] = b[];
    foreach_dimension()
      res[] -= (da.x[1] - da.x[])/Delta;
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  _tol_prev = maxres; //Update for next iteration
#else
  foreach (reduction(max:maxres)) {
    res[] = b[];
    foreach_dimension()
      res[] += (b[1] + b[-1])/10.;
    res[] += (b[1,1] + b[-1,-1] + b[1,-1] + b[-1,1])/100.;
    res[] -= 6./(5.*sq(Delta))*(-4.*a[] +
				4./5.*(a[1] +  a[-1] + a[0,1] + a[0,-1]) +
				1./5.*(a[1,1] + a[1,-1] + a[-1,-1] + a[-1,1]));
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif // !TREE
  
  return maxres;
}
/**
On equidistant grids, the `solve()` function is much faster.
 */
mgstats fast_poisson (scalar a, scalar b) {
  scalar rhs[];
  foreach() {
    rhs[] = b[];
    foreach_dimension() 
      rhs[] += (b[1] + b[-1])/10.; 
    rhs[] += (b[1,1] + b[1,-1] + b[-1,-1] + b[-1,1])/100.;
  }
  solve (a, 6./(5.*sq(Delta))*
	 (-4.   *a[] +
	  4./5.*(a[1] +  a[-1] + a[0,1] + a[0,-1]) +
	  1./5.*(a[1,1] + a[1,-1] + a[-1,-1] + a[-1,1])),
	 rhs);
  return solve_stats;
}

trace
mgstats poisson_compact_mg (struct Poisson p) {
  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;
  scalar a = p.a, b = p.b;
  mgstats s;
#if TREE
  if (tree_is_full())
    s = fast_poisson (a, b);
  else {
    face vector da[]; //Fixme: box boundary conditions needed!
    foreach_dimension()
      da.x.prolongation = refine_face_4_x;
    _tol_prev = -1;
    foreach_face()
      da.x[] = (a[] - a[-1])/Delta;
    s = mg_solve ({a}, {b, da}, residual_compact, relax_compact,
		  &p, p.nrelax, p.res, minlevel = max(1, p.minlevel));
  }
#else
  s = fast_poisson (a, b);
#endif
  if (p.tolerance)
    TOLERANCE = defaultol;
  return s;
}
## test
* [a test](compact_poisson.c)