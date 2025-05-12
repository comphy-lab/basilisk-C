/**
# Compute the frontal area of a fractal tree

A generated fractal tree is projected onto a plane.

![Iteratively adding branches](projecttree2D/c.mp4)

There is an increasing overlap and horizontal orientation of branches.

~~~gnuplot Compare area against total branch length
set logscale xy
set grid
set xlabel 'Branches'
set key box top left
set ylabel 'see key'
plot 'out' u 1:2 t 'Area R^{-2}', '' u 1:3 t 'Length R^{-1}' 
~~~

~~~gnuplot
reset
set logscale xy
set grid
set xlabel 'length R^{-1}'
set ylabel 'area R^{-2}'
set key top left
plot 'out' u 3:2 t 'data', 5*x**0.62 w l t '0.62 scaling'
~~~
*/

#include "grid/cartesian.h"
#include "treegen.h"
#include "utils.h"

Branch * trees;
scalar f[];
face vector fs[];

int main() {
  L0 = 50;
  X0 = -25;
  N = 400;
  init_grid (N);
  srand(0.);
  double angle = -pi/2.; //Rotate tree
  levels = 3;
  trees = tree_skeleton();
  double tl[nb];
  for (int j = 0; j < nb; j++) {
    coord strt = trees[j].start;
    coord end = trees[j].end;
    double bl = 0;
    foreach_dimension(3)
      bl += sq(strt.x - end.x);
    if (j > 0)
      tl[j] = tl[j - 1] + sqrt(bl);
    else
      tl[0] = sqrt(bl);
    trees[j].start.x = cos(angle)*strt.x - sin(angle)*strt.z;
    trees[j].start.z = 0;
    trees[j].end.x = cos(angle)*end.x - sin(angle)*end.z;
    trees[j].end.z = 0;
  }
  int nt = nb;
  
  for (nb = 1; nb <= nt; nb++) {
    tree_interface (trees, f, fs);
    double A = 0;
    foreach() 
      A += dv()*(1 - f[]);
    printf ("%d %g %g\n",nb, A, tl[nb -1]);
    output_ppm (f, file = "c.mp4", opt = "-r 3");
  }
  free (trees);
}