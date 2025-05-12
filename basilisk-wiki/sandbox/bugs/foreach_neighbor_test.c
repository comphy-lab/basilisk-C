/**
# Test case for foreach_neighbor() in multigrid boundaries */

#include "grid/multigrid3D.h"
#include "run.h"

scalar my_scalar[];
vector u[];

int main () 
{
  origin (-.5*L0, -.5*L0, -.5*L0);
  foreach_dimension() 
    periodic(left);
  N = 64;
  run();
}

/** We initialize the Index to -1 everywhere */

event init (i = 0) {
  foreach() 
    my_scalar[] = -1;
  foreach() 
    foreach_dimension() 
      u.x[] = 1.;
}

/** We test what we have inside foreach_neighbor() */

event testvalue (i = 0)
{

  /** Boundary conditions are not automatic when using locate(). */
  
  boundary ({u, my_scalar});

  /** If we use 8 procs, according to the partitioning in Multigrid, 
  we should have (0,0,0) in proc 7 */
  
  Point point = locate(0., 0., 0.);
  if (point.level > -1)
    foreach_neighbor(1)
      printf ("pid() %d (i%d, j%d, k%d), (x%lf, y%lf, z%lf) "
	      "and u.x %lf but my_scalar %d\n", 
	      pid(), point.i, point.j, point.k, x, y, z, 
	      u.x[], (int) my_scalar[]);
}
