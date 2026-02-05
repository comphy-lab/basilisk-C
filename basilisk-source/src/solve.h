/**
# Helper macros to invert (linear) spatial operators

The macros below can be used to easily invert linear systems described
by stencils i.e.
$$
\mathcal{L}(a) = b
$$
For example, let us consider the Poisson equation
$$
\nabla^2 a = b
$$
where $a$ is unknown and $b$ is given. This can be discretised as

~~~literatec
(a[1] + a[-1] + a[0,1] + a[0,-1] - 4.*a[])/sq(Delta) = b[];
~~~

This can be solved using the [`solve()`](#solve) macro below and a
[multigrid solver](poisson.h#mg_solve) with

~~~literatec
solve (a, (a[1] + a[-1] + a[0,1] + a[0,-1] - 4.*a[])/sq(Delta), b[]);
~~~

The macro can take the same optional arguments as
[`mg_solve()`](poisson.h#mg_solve) to tune the multigrid solver.

The macro returns [multigrid statistics](poisson.h#mgstats). 

The [`msolve()`](#msolve) macro generalizes this to systems of linear
equations, with multiple unknown fields. For example let us consider
the coupled reaction-diffusion equations
$$
\begin{aligned}
\partial_tC_1 & = \mu_1\nabla^2C_1 + k_1 C_2, \\
\partial_tC_2 & = \mu_2\nabla^2C_2 + k_2 C_1
\end{aligned}
$$
A first-order, implicit-in-time discretisation could be
$$
\begin{aligned}
\frac{C_1(t) - C_1(t+\Delta t)}{\Delta t} + \mu_1\nabla^2C_1(t+\Delta t) + 
k_1 C_2(t+\Delta t) & = 0, \\
\frac{C_2(t) - C_2(t+\Delta t)}{\Delta t} + \mu_2\nabla^2C_2(t+\Delta t) + 
k_2 C_1(t+\Delta t) & = 0
\end{aligned}
$$
This can be solved using

~~~literatec
scalar C1n[], C2n[];
foreach()
  C1n[] = C1[], C2n[] = C2[];
restriction ({C1n, C2n});

#define laplacian(a) ((a[1,0] + a[-1,0] + a[0,1] + a[0,-1] - 4.*a[])/sq(Delta))

msolve ({C1, C2}, {
  (C1n[] - C1[])/dt + mu1*laplacian(C1) + k1*C2[],
  (C2n[] - C2[])/dt + mu2*laplacian(C2) + k2*C1[]
});
~~~

where the call to `restriction()` ensures that the values of `C1n` and
`C2n` are defined on all the levels of the multigrid. 

Note that the previous example (with only one equation and one
unknown) could also have been solved using

~~~literatec
restriction ({b});
msolve ({a}, {(a[1] + a[-1] + a[0,1] + a[0,-1] - 4.*a[])/sq(Delta) - b[]});
~~~

## Single unknown

Note that the macro below is a slightly simplified version of
the [mg_solve()](poisson.h#mg_solve) and
[mg_cycle()](poisson.h#mg_cycle) functions where more
documentation can be found. */

#include "poisson.h"

macro
mgstats solve (scalar a, double func, double rhs,
	       int nrelax = 4,
	       int minlevel = 0,
	       double tolerance = TOLERANCE)
{{
  mgstats _s = (mgstats){0};
  scalar _res[], _da[];
  scalar_clone (_da, a);
  for (int b = 0; b < nboundary; b++)
    _da.boundary[b] = _da.boundary_homogeneous[b];
  _s.nrelax = nrelax;
  double _resb;
  {
    double maxres = 0.;
    foreach (reduction(max:maxres)) {
      _res[] = rhs - func;
      if (fabs (_res[]) > maxres)
	maxres = fabs (_res[]);
    }
    _resb = _s.resb = _s.resa = maxres;
  }
  for (_s.i = 0; _s.i < NITERMAX && (_s.i < NITERMIN || _s.resa > tolerance); _s.i++) {
    {
      restriction ({_res});
      int _maxlevel = grid->maxdepth;
      int _minlevel = min (minlevel, _maxlevel);
      for (int l = _minlevel; l <= _maxlevel; l++) {
	if (l == _minlevel)
	  foreach_level_or_leaf (l)
	    foreach_blockf (_da)
	      _da[] = 0.;
	else
	  foreach_level (l)
	    foreach_blockf (_da)
	    _da[] = bilinear (point, _da);
	boundary_level ({_da}, l);
	for (int i = 0; i < _s.nrelax; i++) {
	  scalar a = _da;
	  foreach_level_or_leaf (l) {
	    a[] = 0.;
	    double _n = _res[] - func, _d;
	    diagonalize(a)
	      _d = func;
	    a[] = _n/_d;
	  }
/**
Note that the `diagonalize()` operator is not necessary, one could have written instead

~~~literatec
foreach_level_or_leaf (l) {
  a[] = 0.;
  double _d = - func, _n = _res[] + _d;
  a[] = 1.;
  _d += func;
  a[] = _n/_d;
}
~~~
*/
	  boundary_level ({_da}, l);
	}
      }      
      foreach()
	foreach_blockf (a)
	  a[] += _da[];
    }
    {
      double maxres = 0.;
      foreach (reduction(max:maxres)) {
	_res[] = rhs - func;
	if (fabs (_res[]) > maxres)
	  maxres = fabs (_res[]);
      }
      _s.resa = maxres;
    }
    if (_s.resa > tolerance) {
      if (_resb/_s.resa < 1.2 && _s.nrelax < 100)
	_s.nrelax++;
      else if (_resb/_s.resa > 10 && _s.nrelax > 2)
	_s.nrelax--;
    }
    _resb = _s.resa;
  }
  _s.minlevel = minlevel;
  if (_s.resa > tolerance)
    fprintf (stderr,
	     "src/solve.h:%d: warning: convergence for %s not reached after %d iterations\n"
	     "  res: %g nrelax: %d\n", LINENO, a.name,
	     _s.i, _s.resa, _s.nrelax),
      fflush (ferr);
  return _s;
}}

/**
## Multiple unknowns

The `msolve()` macro takes a list of unknowns and a list of equations,
which must have the same length. It also takes the same optional
arguments as [`mg_solve()`](poisson.h#mg_solve) to tune the multigrid
solver. Note however that it does not return the multigrid statistics
but write them in the optional variable given by `mg`. This is
necessary to be able to use the "ellipsis" block which can optionally
be passed to the macro.

This optional ellipsis block is called after each iteration of the
multigrid solver and can be used, for example, to update the
coefficients of the system within the solution procedure (for example
to deal with non-linear couplings).

The overall solution strategy is similar to that used in
[`solve()`](#solve) and [mg_solve()](poisson.h#mg_solve). Aside from
the generalisation to multiple unknowns, the main difference is that
this macro does not require manually splitting the homogeneous and the
non-homogeneous parts of the equations (i.e. the `func` and `rhs`
arguments of [`solve()`](#solve) respectively). */

macro msolve (scalar * X, double * equations,
	      mgstats * mg = NULL,
	      int nrelax = 4,
	      int minlevel = 0,
	      double tolerance = 1e-3)
{{
  mgstats _s = (mgstats){0};
  _s.nrelax = nrelax;
  scalar * _list = (scalar *) X;
  scalar * _lres = list_clone (_list), * _lds = list_clone (_list), * _lrhs = list_clone (_list);
  int _len = list_len (_list);
  double _resb;

  /**
  We compute the initial residuals in `_lres`, for each equation, and
  store the initial solution in `_lds`. The ellipsis block is inserted
  before this evaluation. */
  
  {
    {...}
    double maxres = 0.;
    foreach (reduction(max:maxres)) {
      double * R = equations;
      int i = 0;
      for (scalar res in _lres) {
	res[] = R[i++];
	if (fabs (res[]) > maxres)
	  maxres = fabs (res[]);
      }
      scalar s, ds;
      for (s, ds in _list, _lds)
	ds[] = s[];
    }
    _resb = _s.resb = _s.resa = maxres;
  }

  /**
  On each level of the multigrid hierarchy, we store in `_lrhs` the
  non-homogeneous part of each equation. This is done by resetting the
  vector of unknowns to zero, and re-evaluating each equation. */
  
  int _maxlevel = grid->maxdepth;
  int _minlevel = min (minlevel, _maxlevel);
  reset (_list, 0.);
  for (int l = _minlevel; l <= _maxlevel; l++)
    foreach_level (l) {
      double * R = equations;
      int i = 0;
      for (scalar rhs in _lrhs)
	rhs[] = R[i++];
    }

  /**
  This is the main multigrid iteration loop. */
  
  for (_s.i = 0; _s.i < NITERMAX && (_s.i < NITERMIN || _s.resa > tolerance); _s.i++) {

    /**
    We need homogenous boundary conditions on the vector of unknowns,
    in order to compute the correction to the solution. */
    
    for (scalar s in _list)
      for (int b = 0; b < nboundary; b++)
	s.boundary[b] = s.boundary_homogeneous[b];

    /**
    After having restricted the residual on all levels, we iterate
    from the coarsest to the finest level. On the coarsest level the
    initial guess is zero. On all other levels it is
    bilinearly-interpolated from the coarser level. */
    
    restriction (_lres);
    for (int _l = _minlevel; _l <= _maxlevel; _l++) {
      if (_l == _minlevel)
	foreach_level_or_leaf (_l)
	  for (scalar s in _list)
	    s[] = 0.;
      else
	foreach_level (_l)
	  for (scalar s in _list)
	    s[] = bilinear (point, s);
      boundary_level (_list, _l);

      /**
      This is the relaxation loop. */
      
      for (int _iter = 0; _iter < _s.nrelax; _iter++) {
	foreach_level_or_leaf (_l) {

	  /**
	  We first evaluate the off-diagonal and non-homogeneous terms
	  `_R` of each equation, by setting all diagonal unknowns to
	  zero. */
	  
	  for (scalar _s in _list)
	    _s[] = 0.;
	  double * _R = equations, _D[_len][_len];

	  /**
	  We then compute the matrix `_D` of diagonal coefficients for
	  each unknown and each equation, by setting each diagonal
	  unknown to one for each equation. To get only the diagonal
	  coefficient we need to substract the `_R` vector we just
	  computed. */
	  
	  int _k = 0;
	  for (scalar _s in _list) {
	    _s[] = 1.;
	    double * _r = equations;
	    for (int j = 0; j < _len; j++)
	      _D[_k][j] = _r[j] - _R[j];
	    _s[] = 0.; _k++;
	  }

	  /**
	  The residual for each equation is then computed and stored
	  in `_R`. */
	  
	  _k = 0;
	  scalar rhs, res;
	  for (rhs, res in _lrhs, _lres)
	    _R[_k++] += res[] - rhs[];

	  /**
	  We then solve the resulting system of equations for each
	  diagonal unknown. Note that for the moment we are limited to
	  a maximum of two unknowns, but it should be easy to
	  generalise using e.g. Gauss pivoting etc. */
	  
	  if (_len == 1) {
	    scalar x0 = _list[0];
	    assert (_D[0][0] != 0.);
	    x0[] = - _R[0]/_D[0][0];
	  }
	  else if (_len == 2) {
	    double det = _D[0][0]*_D[1][1] - _D[0][1]*_D[1][0];
	    assert (det != 0.);
	    scalar x0 = _list[0], x1 = _list[1];
	    x0[] = (_D[1][0]*_R[1] - _D[1][1]*_R[0])/det;
	    x1[] = (_D[0][1]*_R[0] - _D[0][0]*_R[1])/det;
	  }
	  else
	    assert (false); // not implemented yet
	}
	boundary_level (_list, _l);
      }
    }

    /**
    Once we have computed the correction on the finest level, we
    update the solution and revert to the original (non-homogeneous)
    boundary conditions. */
    
    foreach() {
      scalar s, ds;
      for (s, ds in _list, _lds)
	s[] += ds[], ds[] = s[];
    }
    {
      scalar s, ds;
      for (s, ds in _list, _lds)
	for (int b = 0; b < nboundary; b++)
	  s.boundary[b] = ds.boundary[b];
    }

    /**
    We then compute the new residual after having applied the optional
    ellipsis block. */
    
    {...}
    double maxres = 0.;
    foreach (reduction(max:maxres)) {
      double * _R = equations;
      int i = 0;
      for (scalar res in _lres) {
	res[] = _R[i++];
	if (fabs (res[]) > maxres)
	  maxres = fabs (res[]);
      }
    }

    /**
    We tune the number of relaxations, based on the convergence rate. */
    
    _s.resa = maxres;
    if (_s.resa > tolerance) {
      if (_resb/_s.resa < 1.2 && _s.nrelax < 100)
	_s.nrelax++;
      else if (_resb/_s.resa > 10 && _s.nrelax > 2)
	_s.nrelax--;
    }
    _resb = _s.resa;
  }

  /**
  We check for convergence, cleanup and return the multigrid statistics. */
  
  _s.minlevel = minlevel;
  if (_s.resa > tolerance) {
    fprintf (stderr, "src/solve.h:%d: warning: convergence for {", LINENO);
    for (scalar a in _list)
      fprintf (stderr, "%s%s", a.name, --_len ? "," : "");
    fprintf (stderr, "} not reached after %d iterations\n"
	     "  res: %g nrelax: %d\n",
	     _s.i, _s.resa, _s.nrelax);
    fflush (stderr);
  }
  delete (_lres), free (_lres);
  delete (_lds), free (_lds);
  delete (_lrhs), free (_lrhs);
  if (((long)mg) + 1 != 1) // just to avoid a warning
    *mg = _s;
}}
