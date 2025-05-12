/**
# Vertex Poisson solver


code adapted from 
http://basilisk.fr/sandbox/Antoonvh/nodal-poisson.h

Special attention is required for dealing with resolution boundaries
where vertices are shared and the restriciton of the residual.
*/
#include "inner-vertex.h"
#include "Antoonvh/my_vertex.h"
#include "poisson.h"

mgstats vpoisson (struct Poisson p) {
  //setup
  vertex scalar da[], res[], a = p.a, b = p.b;
  scalar_clone (da, a);
  da.restriction = restriction_vert;
  da.prolongation = refine_vert;

  vertex scalar mask[];
  mask.restriction = restriction_vert;
  mask.prolongation = refine_vert;

  mask[left]   = 0.;
  mask[right]  = 0.;
  mask[top]    = 0.;
  mask[bottom] = 0.;


  /* da[left]   = 0.; */
  /* da[right]  = 0.; */
  /* da[top]    = 0.; */
  /* da[bottom] = 0.; */

  foreach_vertex()
    mask[] = 1.;
  boundary ({mask});
  restriction({mask});
  for (int l = 0; l <= depth(); l++) {
    boundary_level({mask}, l);

  }

  if (p.res)
    res = p.res[0];
  else 
    scalar_clone (res, b);
//  res.restriction = restriction_vert;
  res.prolongation = refine_vert;

  // BD: In principle, there is no need to define a BC for the rhs (or
  // residual).  However, for the Multigrid method, these BC are necessary to
  // converge in case a and b do not have the same BC (e.g. for the QG model with no slip BC)
  res[left]   = 0.;
  res[right]  = 0.;
  res[top]    = 0.;
  res[bottom] = 0.;


  mgstats mg; mg.sum = HUGE, mg.resa = HUGE;
  mg.nrelax = p.nrelax ? p.nrelax : 5;
  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;
  mg.minlevel = p.minlevel ? p.minlevel : 1;
  /* Solver Iterations */
  for (mg.i = 0; mg.i < NITERMAX; mg.i++) {
    // Residual
    double max = 0;
    foreach_vertex(reduction (max:max)) {
      res[] = b[]*mask[];
      foreach_dimension() {
          res[] -= (a[-1] - 2.*a[] + a[1])/(sq(Delta))*mask[];
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
    for (int l = 0; l <= depth(); l++) {
      boundary_level({res}, l);      
    }

    // Guess
    foreach_vertex_level(mg.minlevel)
      da[] = 0;
    // Up-cycle
    for (int l = mg.minlevel; l <= depth(); l++) {
      boundary_level({da}, l);
      // Relaxation sweep
      for (int rel = 0; rel < mg.nrelax; rel++) {
        foreach_vertex_level(l) {
          double d = 0;
          da[] = -res[]*sq(Delta);
          foreach_dimension() {
              da[] += (da[1] + da[-1])*mask[];
              d += 2.;
          }
          da[] /= d;
        }
        boundary_level({da}, l);
      }
      // Prolongation
      //    printf("LEVEL: %d\n",l);

      
      // modif BD
      if (l < depth()){
        foreach_vertex_level(l){
          //   printf("%g \t %g \t %g\n",x,y,da[]);
          refine_vert (point, da);
        }
        boundary_level({da}, l+1);

      }
    }
    // Correction
    foreach_vertex()
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
