/**
# Vertex Poisson solver


code adapted from 
http://basilisk.fr/sandbox/Antoonvh/nodal-poisson.h

Special attention is required for dealing with resolution boundaries
where vertices are shared and the restriciton of the residual.
*/
#include "inner-vertex_mpi.h"
#include "my_vertex.h"
#include "poisson.h"

mgstats vpoisson (struct Poisson p) {
  //setup
  vertex scalar da[], res[], a = p.a, b = p.b;
  scalar_clone (da, a);
  da.restriction = restriction_vert;
  da.prolongation = refine_vert;

  if (p.res)
    res = p.res[0];
  else 
    scalar_clone (res, b);
//  res.restriction = restriction_vert;
  res.prolongation = refine_vert;
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
    foreach_inner_vertex(reduction (max:max)) {
      res[] = b[];
      foreach_dimension() {
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
    foreach_vertex_level(mg.minlevel)
      da[] = 0;
    // Up-cycle
    for (int l = mg.minlevel; l <= depth(); l++) {
      boundary_level({da}, l);
      // Relaxation sweep
      for (int rel = 0; rel < mg.nrelax; rel++) {
        foreach_inner_vertex_level(l) {
          double d = 0;
          da[] = -res[]*sq(Delta);
          foreach_dimension() {
              da[] += (da[1] + da[-1]);
              d += 2.;
          }
          da[] /= d;
        }
        boundary_level({da}, l);
      }
      // Prolongation
      //    printf("LEVEL: %d\n",l);

      
      // modif BD
      if (l < depth())
        foreach_inner_vertex_level(l){
          //   printf("%g \t %g \t %g\n",x,y,da[]);
          refine_vert (point, da);
        }
    }
    // Correction
    foreach_inner_vertex()
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
