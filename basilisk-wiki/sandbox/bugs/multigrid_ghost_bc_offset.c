//#include "grid/quadtree.h"
#include "grid/multigrid.h"

/**
When using 2 layers of Ghost cells in the multigrid solver, the stencil used when evaluating the value of the ghost cell is offset away from the boundary at the second ghost cell. (Did not test for more than 2 ghost cells)

1rst Ghost cell -> centered at cell adjacent to boundary

2nd Ghost cell -> centered at 2nd cell from boundary

This occurs when using the multigrid solver but not the tree solver. 

Note that the same function defined by s[boundary] is used to evaluate all layers of ghost cells.

*/

#define BGHOSTS 2

int main() {
  
  init_grid (4);
  L0=4;

  scalar s[];
  foreach()
    s[]=x+0.5;
  s[right]=s[];
  s[left]=s[];
  boundary ((scalar *){s});

  foreach_boundary (left){
    fprintf (stdout, "%g %g %g %g\n",s[-2],s[-1],s[],s[1]);
  }
  foreach_boundary (right){
    fprintf (stdout, "%g %g %g %g\n",s[2],s[1],s[],s[-1]);
  }
}

  /**
s is 1,2,3,4 horizonally and uniform vertically

s[left ] should be 1,1 but is instead 2,1

s[right] should be 4,4 but is instead 4,3
   */