/**
# Compute the boundary conditions at second order */

#include "grid/multigrid.h"

#define BGHOSTS 2

@define n_x neighbor.i
@define n_y neighbor.j
@define n_z neighbor.k
foreach_dimension() {
  @define layer_nr_x (n_x < GHOSTS ? (GHOSTS - n_x) : -(n_x - (1 << level) - GHOSTS + 1)) 
  @define bc_x(a)    ((a) - sign(layer_nr_x)*(-Delta/2. + fabs (layer_nr_x)*Delta))
}

int main()
{
  L0 = 4;
  N = 4;
  init_grid (N);
  output_cells (stdout);

  scalar s[];
  s[left] = (bc_x(x)); // At the boundary x=0, therefore the cell center is x - Delta/2.
  s[right] = 2.*(bc_x(x));
  s[bottom] = 2.*(bc_y(y));
  s[top] = (bc_y(y));
  foreach()
    s[] = 0.;
  boundary ((scalar *) {s});
  
  foreach()
    fprintf (stderr, "%g %g %g\n", x, y, s[]);

  foreach_boundary(left)
    for (int i = -2; i < 0; i++)
      fprintf (stderr, "%g %g %g\n", x + i*Delta, y, s[i]);
  foreach_boundary(right)
    for (int i = 1; i <= 2; i++)
      fprintf (stderr, "%g %g %g\n", x + i*Delta, y, s[i]);

  foreach_boundary(bottom)
    for (int i = -2; i < 0; i++)
      fprintf (stderr, "%g %g %g\n", x, y + i*Delta, s[0,i]);
  foreach_boundary(top)
    for (int i = 1; i <= 2; i++)
      fprintf (stderr, "%g %g %g\n", x, y + i*Delta, s[0,i]);
}

/**
# Results

~~~gnuplot
unset key 
unset border
unset tics
unset xlabel
unset ylabel
set size ratio -1
plot 'out' w l lc rgb "black", \
     'log' u 1:2:3 with labels tc rgb "blue"
~~~
*/
