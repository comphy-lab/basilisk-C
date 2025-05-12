/**
# Box-boundary values are not set consistently

Here the quadtree version is demonstrated. The bottom-left ghost cell
has a value of `s = 5`, However, on a multigrid (or Cartesian grid) it becomes `s = 1`. I wish it was 3?
 
~~~gnuplot Values and their locations
set xr [-0.3:0.6]
set yr [-0.3:0.6]
set size square
set key off
plot 'log' w l lw 2, 'out' u 1:2:3 w labels
~~~
 */
//#include "grid/multigrid.h" // <- uncomment for MG

scalar s[];

s[left]   = dirichlet (2.);
s[bottom] = dirichlet (2.);

int main() {
  init_grid (4);
  foreach()
    s[] = 1.;
  boundary ({s});
  output_cells (stderr);
  foreach() {
    foreach_neighbor(1)
      printf ("%g %g %g\n", x, y, s[]);
    return 0;
  }
}
