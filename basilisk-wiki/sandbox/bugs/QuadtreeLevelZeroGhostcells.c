/**
# Inconsistent 5x5 stencil on root level of trees */

#define BGHOSTS 2

int main() {
  init_grid(2);

  scalar A[];
  foreach_level(0)
    A[] = 1.;

  foreach_dimension() {
    A[left] = neumann (0);
    A[right] = neumann (0);
  }
  boundary_level({A}, 0);

  foreach_level(0)
    foreach_neighbor()
      fprintf (stderr, "%g %g %g\n", x, y, A[]);
}

/**
There are NaNs in some locations.

~~~gnuplot
set size ratio -1
plot 'log' u 1:2:3 w labels
~~~
*/
