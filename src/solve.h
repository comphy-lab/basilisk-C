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

The [multigrid statistics](poisson.h#mgstats) are stored in
`solve_stats`. */

#include "poisson.h"
static mgstats solve_stats;

/**
## Implementation

Note that the large macro below is a slightly simplified version of
the [mg_solve()](poisson.h#mg_solve) and
[mg_cycle()](poisson.h#mg_cycle) functions where more
documentation can be found. */

#define solve(_a, _func, _rhs, ...)					\
  do {									\
    solve_stats = (mgstats){0};						\
    struct { int dummy, nrelax, minlevel; double tolerance; }		\
      p = {0, __VA_ARGS__};						\
    scalar _res[], _da[];						\
    scalar_clone (_da, _a);						\
    for (int b = 0; b < nboundary; b++)					\
      _da.boundary[b] = _da.boundary_homogeneous[b];			\
    solve_stats.nrelax = p.nrelax > 0 ? p.nrelax : 4;			\
    double resb;							\
    {									\
      double maxres = 0.;						\
      foreach (reduction(max:maxres)) {					\
	_res[] = _rhs[] - (_func);					\
	if (fabs (_res[]) > maxres)					\
	  maxres = fabs (_res[]);					\
      }									\
      resb = solve_stats.resb = solve_stats.resa = maxres;		\
    }									\
    if (p.tolerance == 0.)						\
      p.tolerance = TOLERANCE;						\
    for (solve_stats.i = 0;						\
	 solve_stats.i < NITERMAX &&					\
	   (solve_stats.i < NITERMIN || solve_stats.resa > p.tolerance); \
	 solve_stats.i++) {						\
      {									\
	restriction ({_res});						\
	int maxlevel = grid->maxdepth;					\
	int minlevel = min (p.minlevel, maxlevel);			\
	for (int l = minlevel; l <= maxlevel; l++) {			\
	  if (l == minlevel)						\
	    foreach_level_or_leaf (l)					\
	      foreach_blockf (_da)					\
	        _da[] = 0.;						\
	  else								\
	    foreach_level (l)						\
		foreach_blockf (_da)					\
		  _da[] = bilinear (point, _da);			\
	  boundary_level ({_da}, l);					\
	  for (int i = 0; i < solve_stats.nrelax; i++) {		\
	    scalar _a = _da;						\
	    foreach_level_or_leaf (l) {					\
	      _a[] = 0.;						\
	      double n = _res[] - (_func), d;				\
	      diagonalize(_a) {						\
		d = (_func);						\
	      }								\
	      _a[] = n/d;						\
	    }								\
	    boundary_level ({_da}, l);					\
	  }								\
	}								\
	foreach()							\
	  foreach_blockf (_a)						\
	    _a[] += _da[];						\
      }									\
      {									\
	double maxres = 0.;						\
	foreach (reduction(max:maxres)) {				\
	  _res[] = _rhs[] - (_func);					\
	  if (fabs (_res[]) > maxres)					\
	    maxres = fabs (_res[]);					\
	}								\
	solve_stats.resa = maxres;					\
      }									\
      if (solve_stats.resa > p.tolerance) {				\
	if (resb/solve_stats.resa < 1.2 && solve_stats.nrelax < 100)	\
	  solve_stats.nrelax++;						\
	else if (resb/solve_stats.resa > 10 && solve_stats.nrelax > 2)	\
	  solve_stats.nrelax--;						\
      }									\
      resb = solve_stats.resa;						\
    }									\
    solve_stats.minlevel = p.minlevel;					\
    if (solve_stats.resa > p.tolerance)					\
      fprintf (ferr,							\
	       "WARNING: convergence for %s not reached "		\
	       "after %d iterations\n"					\
	       "  res: %g nrelax: %d\n", _a.name,			\
	       solve_stats.i, solve_stats.resa, solve_stats.nrelax),	\
	fflush (ferr);							\
  } while (0)
