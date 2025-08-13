/**
# Solve multiple Poisson problems in the same cycle.

This is a copy-paste of the regular Poisson solver, adjusted to solve multiple poisson problems defined by lists of `a`, `b`, (same Alpha, Lambda)
 */
#include "poisson.h"

static void relax_list (scalar * al, scalar * bl, int l, void * data)
{
  struct Poisson * p = (struct Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  
  /**
  We use either Jacobi (under)relaxation, Gauss-Seidel or we directly
  reuse values as soon as they are updated. For Jacobi, we need to
  allocate space for the new field *c*. Jacobi is useful mostly as it
  gives results which are independent of the order in which the cells
  are traversed. This is not the case for the simple traversal, which
  means for example that results will depend on whether a tree or a
  multigrid is used (because cells will be traversed in a different
  order). The same comment applies to OpenMP or MPI parallelism. In
  practice however Jacobi convergence tends to be slower than simple
  reuse. */
  scalar a, b, c;
#if JACOBI
  scalar * cl = list_clone(al);

#else
  scalar * cl = al;
#endif
  /**
  On GPUs, we use red/black Gauss-Seidel relaxation, which requires
  two loops (for odd/even indices). Note also that, unlike the other
  option, red/black relaxation should be deterministic. */
  
#if GAUSS_SEIDEL || _GPU
  for (int parity = 0; parity < 2; parity++)
    foreach_level_or_leaf (l, nowarning)
      if (level == 0 || ((point.i + parity) % 2) != (point.j % 2))
#else
  foreach_level_or_leaf (l, nowarning)
#endif
  {

    /**
    We use the face values of $\alpha$ to weight the gradients of the
    5-points Laplacian operator. We get the relaxation function. */

    for (a, b, c in al, bl, cl) {
      double n = - sq(Delta)*b[], d = - lambda[]*sq(Delta);
      foreach_dimension() {
	n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
	d += alpha.x[1] + alpha.x[];
      }
#if EMBED
      if (p->embed_flux) {
	double c, e = p->embed_flux (point, a, alpha, &c);
	n -= c*sq(Delta);
	d += e*sq(Delta);
      }
      if (!d)
	c[] = 0., b[] = 0.;
      else
#endif // EMBED
	c[] = n/d;
    }
  }

  /**
     For weighted Jacobi we under-relax with a weight of 2/3. */
  
#if JACOBI
  foreach_level_or_leaf (l) {
    for (scalar a, c in al, cl)  
      a[] = (a[] + 2.*c[])/3.;
  }
#endif
  
#if TRASH
  for (scalar a in al) {
    scalar a1[];
    foreach_level_or_leaf (l)
      a1[] = a[];
    trash ({a});
    foreach_level_or_leaf (l)
      a[] = a1[];
  }
#endif
#if JACOBI
  delete (cl);
  free (cl);
#endif
}


static double residual_list (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a, b, res;
  face vector g;
  struct Poisson * p = (struct Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector * gl = NULL;
  for (scalar a in al) {
    face vector g = new_face_vector("g");
    gl = vectors_append (gl, g);
  }
  foreach_face() {
    for (a, g in al, gl)      
      g.x[] = alpha.x[]*face_gradient_x (a, 0);
  }
  foreach (reduction(max:maxres), nowarning) {
    for (a, b, res, g in al, bl, resl, gl) {
      res[] = b[] - lambda[]*a[];
      foreach_dimension()
	res[] -= (g.x[1] - g.x[])/Delta;
#if EMBED
      if (p->embed_flux) {
	double c, e = p->embed_flux (point, a, alpha, &c);
	res[] += c - e*a[];
      }
#endif // EMBED    
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
    }
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres), nowarning) {
    for (a, b, res, g in al, bl, resl, gl) {
      res[] = b[] - lambda[]*a[];
      foreach_dimension()
	res[] += (alpha.x[0]*face_gradient_x (a, 0) -
		  alpha.x[1]*face_gradient_x (a, 1))/Delta;  
#if EMBED
      if (p->embed_flux) {
	double c, e = p->embed_flux (point, a, alpha, &c);
	res[] += c - e*a[];
    }
#endif // EMBED
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
    }
  }
#endif // !TREE
  
  delete ((scalar*)gl);
  free (gl);
  
  return maxres;
}

trace
mgstats poisson_list (scalar * al, scalar * bl,
		      (const) face vector alpha = {{-1}},
		      (const) scalar lambda = {-1},
		      double tolerance = 0.,
		      int nrelax = 4,
		      int minlevel = 0,
		      scalar * res = NULL,
		      double (* flux) (Point, scalar, vector, double *) = NULL)
{
  /**
  If $\alpha$ or $\lambda$ are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (alpha.x.i < 0)
    alpha[] = {1.,1.,1.};
  if (lambda.i < 0)
    lambda[] = 0.;
  
  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  restriction ({alpha,lambda});

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (tolerance)
    TOLERANCE = tolerance;

  struct Poisson p = {al[0], bl[0], alpha, lambda, tolerance, nrelax, minlevel, res };
#if EMBED
  if (!flux && a.boundary[embed] != symmetry)
    p.embed_flux = embed_flux;
  else
    p.embed_flux = flux;
#endif // EMBED
  mgstats s = mg_solve (al, bl, residual_list, relax_list, &p,
			nrelax, res, max(1, minlevel));

  /**
  We restore the default. */

  if (tolerance)
    TOLERANCE = defaultol;

  return s;
}


