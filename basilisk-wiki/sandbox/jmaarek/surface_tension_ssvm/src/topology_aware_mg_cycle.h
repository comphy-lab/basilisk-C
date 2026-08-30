/**
# topology_aware_mg_cycle.h

Simple impmentation of a topologically aware multigrid. 
Relative to the default multigrid cycle the only thing that changes is that
if we are in a near interface cell we use direct injection as opposed to 
bilinear prolongation. For poisson problems with large density ratios, this is observed to speed up convergence potentially
signficantly. This can be explained by direct injection preventing blending 
across a discontinuity.

I tried implementing  the same for the viscous term but found no signficant
speed up. In general it was slightly worse in this case than the standard approach
for the case i tried. 

Defines an independently-named parallel chain -- mg_cycle_topology_aware,
mg_solve_topology_aware, poisson_topology_aware, project_topology_aware,
-- none of which share a name with anything in
poisson.h, so there is no possibility of a redefinition
error the way there was trying to override mg_cycle/mg_solve directly
(see the discussion earlier in this conversation for why that approach
can't work against an unconditionally-defined function).



*/

#include "poisson.h"


extern scalar f;

static inline double bilinear_topology_aware (Point point, scalar s)
{
  double f0 = coarse (f);
  double fmin = f0, fmax = f0;
 
  double f1 = coarse (f, child.x);
  if (f1 < fmin) fmin = f1;
  if (f1 > fmax) fmax = f1;
 
#if dimension >= 2
  double f2 = coarse (f, 0, child.y);
  if (f2 < fmin) fmin = f2;
  if (f2 > fmax) fmax = f2;
  double f3 = coarse (f, child.x, child.y);
  if (f3 < fmin) fmin = f3;
  if (f3 > fmax) fmax = f3;
#endif
#if dimension >= 3
  double f4 = coarse (f, 0, 0, child.z);
  if (f4 < fmin) fmin = f4;
  if (f4 > fmax) fmax = f4;
  double f5 = coarse (f, child.x, 0, child.z);
  if (f5 < fmin) fmin = f5;
  if (f5 > fmax) fmax = f5;
  double f6 = coarse (f, 0, child.y, child.z);
  if (f6 < fmin) fmin = f6;
  if (f6 > fmax) fmax = f6;
  double f7 = coarse (f, child.x, child.y, child.z);
  if (f7 < fmin) fmin = f7;
  if (f7 > fmax) fmax = f7;
#endif
 
  if (fmin > 1. - 1e-6 ||   // all pure liquid: even the smallest
			     // corner value is (numerically) 1
      fmax < 1e-6)          // all pure gas: even the largest
			     // corner value is (numerically) 0
    return bilinear (point, s);   // single phase throughout this
				   // coarse neighbourhood -- no
				   // interface anywhere inside it, safe
				   // to blend
  return coarse (s);              // interface present somewhere in
				   // this coarse neighbourhood --
				   // fall back to injection
}



trace
void mg_cycle_topology_aware (scalar * a, scalar * res, scalar * da,
			       void (* relax) (scalar * da, scalar * res,
					       int depth, void * data),
			       void * data,
			       int nrelax, int minlevel, int maxlevel)
{

  restriction (res);

  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {

    if (l == minlevel)
      foreach_level_or_leaf (l)
	for (scalar s in da)
	  foreach_blockf (s)
	    s[] = 0.;
    else {
      boundary_level (da, l - 1);
      foreach_level (l)
	for (scalar s in da)
	  foreach_blockf (s)
	    s[] = bilinear_topology_aware(point, s);   // swap with topology-aware prolongation
    }

    for (int i = 0; i < nrelax; i++) {
      boundary_level (da, l);
      relax (da, res, l, data);
    }
  }

  foreach() {
    scalar s, ds;
    for (s, ds in a, da)
      foreach_blockf (s)
	s[] += ds[];
  }
}

/**
## mg_solve_topology_aware()

Direct copy of poisson.h's own mg_solve, with the one call to mg_cycle
literally written as mg_cycle_topology_aware -- no macro needed here,
since this is a plain function definition under a name that doesn't
exist anywhere else in the program.
*/

trace
mgstats mg_solve_topology_aware (scalar * a, scalar * b,
				  double (* residual) (scalar * a, scalar * b,
						       scalar * res, void * data),
				  void (* relax) (scalar * da, scalar * res,
						  int depth, void * data),
				  void * data = NULL,
				  int nrelax = 4,
				  scalar * res = NULL,
				  int minlevel = 0,
				  double tolerance = TOLERANCE)
{
  restriction({f});
  scalar * da = list_clone (a), * pres = res;
  if (!res)
    res = list_clone (b);

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];

  mgstats s = {0};
  double sum = 0.;
  scalar rhs = b[0];
  foreach (reduction(+:sum))
    sum += rhs[];
  s.sum = sum;
  s.nrelax = nrelax > 0 ? nrelax : 4;

  double resb;
  resb = s.resb = s.resa = (* residual) (a, b, res, data);

  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > tolerance);
       s.i++) {
    mg_cycle_topology_aware (a, res, da, relax, data,
			      s.nrelax,
			      minlevel,
			      grid->maxdepth);
    s.resa = (* residual) (a, b, res, data);

    if (s.resa > tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
	s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
	s.nrelax--;
    }

    resb = s.resa;
  }
  s.minlevel = minlevel;

  if (s.resa > tolerance) {
    scalar v = a[0];
    fprintf (ferr,
	     "topology_aware_mg_cycle.h:%d: warning: convergence for %s "
	     "not reached after %d iterations\n"
	     "  res: %g sum: %g nrelax: %d tolerance: %g\n", LINENO, v.name,
	     s.i, s.resa, s.sum, s.nrelax, tolerance), fflush (ferr);
  }

  if (!pres)
    delete (res), free (res);
  delete (da), free (da);

  return s;
}

/**
## poisson_topology_aware()

Direct copy of poisson.h's own poisson(), calling
mg_solve_topology_aware instead of mg_solve. relax/residual referenced
here are poisson.h's own (static, but directly callable from this
header since Basilisk flattens the whole program into one translation
unit -- the same reasoning established earlier in this conversation).
EMBED handling kept for structural fidelity with poisson.h; untested
specifically with EMBED enabled.
*/

trace
mgstats poisson_topology_aware (scalar a, scalar b,
				 (const) face vector alpha = {{-1}},
				 (const) scalar lambda = {-1},
				 double tolerance = 0.,
				 int nrelax = 4,
				 int minlevel = 0,
				 scalar * res = NULL,
				 double (* flux) (Point, scalar, vector, double *) = NULL)
{
  if (alpha.x.i < 0)
    alpha[] = {1., 1., 1.};
  if (lambda.i < 0)
    lambda[] = 0.;

  restriction ({alpha, lambda});

  double defaultol = TOLERANCE;
  if (tolerance)
    TOLERANCE = tolerance;

  struct Poisson p = {a, b, alpha, lambda, tolerance, nrelax, minlevel, res};
#if EMBED
  if (!flux && a.boundary[embed] != symmetry)
    p.embed_flux = embed_flux;
  else
    p.embed_flux = flux;
#endif
  mgstats s = mg_solve_topology_aware ({a}, {b}, residual, relax, &p,
					nrelax, res, max (1, minlevel));

  if (tolerance)
    TOLERANCE = defaultol;

  return s;
}

/**
## project_topology_aware()

Direct copy of poisson.h's own project(), calling poisson_topology_aware
instead of poisson(). This is the function you actually redirect
centered.h's internal call to -- see the wiring example below.
*/

trace
mgstats project_topology_aware (face vector uf, scalar p,
				 (const) face vector alpha = unityf,
				 double dt = 1.,
				 int nrelax = 4)
{
  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
  }

  mgstats mgp = poisson_topology_aware (p, div, alpha,
					 tolerance = TOLERANCE/sq(dt),
					 nrelax = nrelax);

  foreach_face()
    uf.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}

#define project(...) project_topology_aware(__VA_ARGS__)


/**
## Wiring example (goes in your main .c file, NOT in this header)

#include "grid/quadtree.h"
#include "embed.h"
#include "topology_aware_mg_cycle.h"   // this file -- includes
                                        // poisson.h
                                        // for real, first

#include "navier-stokes/centered.h"    // its internal calls to
                                        // project(...) inside its
                                        // own events
#include "two-phase.h"

*/
