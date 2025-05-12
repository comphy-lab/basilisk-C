/**
# Prolongate fluxes 

We check if we can prolongate halo fluxes rather than restrict them.
*/

#include "grid/quadtree.h"
#include "prolongate_halo_flux.h"

face vector f[];
f.t[bottom] = dirichlet (level + 0.1);

int main() {
  init_grid (16);
  unrefine(x < 0.7 && y < 0.7);
  foreach_cell() {
    f.x[1] = f.x[] = level;
    f.y[0,1] = f.y[] = -level;
  }
  /**
Since `boundary()` is too eager to restrict the face values, we need to set the box-boundary-ghost values via this method (:S):
   */
  Boundary * b = qcalloc (1, Boundary);
  b->level = box_boundary_level;
  b->restriction = masked_boundary_restriction;
  add_boundary (b);
  for (int l = 0 ; l <= depth() ; l++)
    box_boundary_level (b, {f.x}, l);
  /**
     Now we may prolongate:
  */
  halo_flux_prolongate({f});
  /**
We check if this has worked as inteded
   */
  FILE * fp = fopen ("faces","w");
  foreach_face()
    fprintf(fp, "%g %g %g %d\n",x, y, f.x[], level);
}
/**
~~~gnuplot If you look closely, everything worked as intented. 
set xr [-0.1:1.1]
set yr [-0.1:1.1]
plot 'faces' u 1:2:3 with labels
~~~
*/