/**
## Staggered Helmholtz-Hodge decomposition in three dimensions.

Extending the Helmholtz-Hodge decomposition example of [Jose-Maria
Fullana](/sandbox/jmf/Hodge/hodge.c), we decompose a face-averaged
vector field ($\vec{F}$) as the sum of a discrete gradient of a
cell-averaged field ($\nabla_D a$) and the exact face-averaged curl of
an edge-averaged vector field ($\nabla \times \vec{A}$). The ansantz
is that the discrete field is defined by its exact cell-averaged
divergence and approximate edge-averaged curl.

~~~gnuplot The reconstruction has been accomplished
set size ratio 1
set xlabel 'Original'
set ylabel 'reconstructed'
set grid
set key top left
plot x t'1:1', 'out' t 'data' 
~~~


~~~gnuplot The reconstruction error is controlled by the TOLERANCE
set size ratio 1
set xlabel 'Value'
set ylabel 'Error'
set grid
set key right outside
set yr [-1.4e-5: 1.4e-5]
plot 'out' u 1:3 t 'data', 1e-5 lw 3 t 'TOLERANCE', -1e-5 lw 3 t 'TOLERANCE'
~~~
*/
#include "grid/octree.h"
#include "utils.h"
#include "poisson.h"
#include "vector_pot.h"
/**
## Test the Vector-potential solver

 */
int main() {
  /**
The case concerns a periodic field of noise with zero means.
   */
  X0 = Y0 = Z0 = -L0/2;
  foreach_dimension()
    periodic(left);
  init_grid (8);
  TOLERANCE = 1e-5;
  face vector F[];
  vector A[]; // Edge-averaged vector field components
  scalar a[], div[];   // cell-averaged fields
  foreach_face()
    F.x[] = noise();
  stats f_x = statsf(F.x);
  stats f_y = statsf(F.y);
  stats f_z = statsf(F.z);
  foreach_face()
    F.x[] -= f_x.sum/(cube(L0));
  /**
Fields are initialized. Among others, the exact cell-averaged
divergence.
   */
  foreach() {
    a[] = 0;
    div[] = 0;
    foreach_dimension() {
      A.x[] = 0;
      div[] += (F.x[1] - F.x[])/Delta;
    }
  }
/**
We solve the for scalar and vector potential.
 */
  poisson (a, div);
  vector_potential (F, A, Coulomb_gauge = true);
  /**
We can print the reconstructed field values and compare against the orignal values.
   */
  foreach_face() {
    double reconst = (A.y[] + A.z[0,1] - A.y[0,0,1] - A.z[] + a[] - a[-1])/Delta;
    printf ("%g %g %g\n", F.x[], reconst , F.x[] - reconst);
  }
}