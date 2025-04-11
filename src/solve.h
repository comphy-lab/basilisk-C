/**
# Helper macro to invert (linear) spatial operators

The macro below can be used to easily invert linear systems described
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

This can be solved using the macro below and a [multigrid
solver](poisson.h#mg_solve) with

~~~literatec
solve (a, (a[1] + a[-1] + a[0,1] + a[0,-1] - 4.*a[])/sq(Delta), b);
~~~

The macro can take the same optional arguments as
[mg_solve()](poisson.h#mg_solve) to tune the multigrid solver.

The macro returns [multigrid statistics](poisson.h#mgstats). */

#include "poisson.h"

/**
## Implementation

Note that the macro below is a slightly simplified version of
the [mg_solve()](poisson.h#mg_solve) and
[mg_cycle()](poisson.h#mg_cycle) functions where more
documentation can be found. */

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
