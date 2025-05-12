/**
# Vertex and face boundary conditions for layered variables

Same example as src/test/vertices-bc.c but with layered variables

*/
#define LAYERS 1

#include "grid/multigrid.h"

vertex scalar omega;
face vector f;

int main()
{
  nl = 3;
  N = 4;
  init_grid (N);

  omega = new vertex scalar[nl];
  f = new face vector[nl];


  omega[left] = 0.;
  omega[right] = 0.;
  omega[top] = 0.;
  omega[bottom] = 0.;

  foreach_vertex() {
    int l = 1;
    foreach_layer()
      omega[] = l++;
  }

  f.n[left] = 0.;
  f.n[right] = 0.;
  f.n[top] = 0.;
  f.n[bottom] = 0.;
  
  foreach_face() {
    int l = 1;
    foreach_layer()
      f.x[] = l++;
  }

  /**
  <div class="warning">
  For the moment, stencils for boundary conditions on vertex fields
  are not automatic, so that the manual call to `boundary()` below is
  necessary. This should be fixed.</div> */
  
  boundary ({omega});
  
  foreach_vertex()
    foreach_layer()
      fprintf (qerr, "v %g\t%g\t%g\n", x, y, omega[]);

  foreach_face()
    foreach_layer()
      fprintf (qerr, "f %g\t%g\t%g\n", x, y, f.x[]);

  delete ({omega, f});


}

/**
~~~gnuplot
unset key
set size ratio -1
set xrange [-0.1:1.1]
set yrange [-0.1:1.1]
set xtics 0,0.25,1
set ytics 0,0.25,1
set grid
plot '< grep f log' u 2:3:4 w labels, \
     '< grep v log' u 2:3:4 w labels
~~~
*/
