/**
# Vertex Poisson solver

Special attention is required for dealing with resolution boundaries
where vertices are shared and the restriciton of the residual.
*/
#include "my_vertex.h"
#include "poisson.h"

mgstats vpoisson (struct Poisson p) {
  //setup
  scalar da[], res[], a = p.a, b = p.b;
  scalar_clone (da, a);
  if (p.res)
    res = p.res[0];
  else 
    scalar_clone (res, b);
  res.prolongation = refine_vert;
  mgstats mg; mg.sum = HUGE, mg.resa = HUGE;
  mg.nrelax = p.nrelax ? p.nrelax : 5;
  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;
  mg.minlevel = p.minlevel ? p.minlevel : 1;
  // Solver Iterations
  for (mg.i = 0; mg.i < NITERMAX; mg.i++) {
    // Residual
    double max = 0;
    foreach(reduction (max:max)) {
      res[] = b[];
      foreach_dimension() {
	if (is_prolongation(neighbor(-1))) // Coarse stencil
	  res[] -= (a[-2] - 2.*a[] + a[2])/(4*sq(Delta));
	else
	  res[] -= (a[-1] - 2.*a[] + a[1])/(sq(Delta));
      }
      if (fabs(res[]) > max)
	max = fabs(res[]);
    }
    // Statistics 
    mg.resa = max;
    if (mg.i == 0) 
      mg.resb = max;
    // Break out
    if (max < TOLERANCE && mg.i >= NITERMIN)  
      break;
    // Residual on levels requires attention
    res.restriction = restriction_vert;
    boundary ({res});
    res.restriction = restriction_coarsen_vert; 
    multigrid_restriction ({res});
    // Guess
    foreach_level(mg.minlevel)
      da[] = 0;
    // Up-cycle
    for (int l = mg.minlevel; l <= depth(); l++) {
      boundary_level({da}, l);
      // Relaxation sweep
      for (int rel = 0; rel < mg.nrelax; rel++) {
	foreach_level_or_leaf(l) {
	  double d = 0;
	  da[] = -res[]*sq(Delta);
	  foreach_dimension() {
	    if (is_prolongation(neighbor(-1))) { // Coarse stencil
	      da[] += (da[2] + da[-2])/4.;
	      d += 0.5;
	    } else {
	      da[] += (da[1] + da[-1]);
	      d += 2.;
	    }
	  }
	  da[] /= d;
	}
	boundary_level({da}, l);
      }
      // Prolongation
      foreach_coarse_level(l) 
	refine_vert (point, da);
    }
    // Correction
    foreach()
      a[] += da[];
    boundary ({a});
  }
  if (mg.resa > TOLERANCE)
    fprintf (stderr, "Convergence for %s not reached.\n"
	     "mg.i = %d, mg.resb: %g mg.resa: %g\n",
	     a.name, mg.i, mg.resb,  mg.resa);
  if (p.tolerance)
    TOLERANCE = defaultol;
  return mg;
}
/**
## Application to the solenoidal projection of a vector field
 */
mgstats vproject (vector v, scalar p, double tol) {
  scalar div[];
  foreach() {
    div[] = 0;
    foreach_dimension()
      div[] += (v.x[1] - v.x[-1])/(2*Delta);
  }
  mgstats mg = vpoisson (p, div, tolerance = tol);
  foreach() {
    foreach_dimension() {
      if (is_prolongation(neighbor(-1))) 
	v.x[] -= (p[2] - p[-2])/(4*Delta);
      else
	v.x[] -= (p[1] - p[-1])/(2*Delta);
    }
  }
  boundary ((scalar*){v});
  return mg;
}

/**
## Tests 

*[Test for the Poisson-solver function](tvpoisson.c)  
*[Test for the projection function](tvproject.c)  

 */
