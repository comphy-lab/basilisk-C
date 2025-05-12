/**
# Solve $\nabla^2 a = b$ for $\nabla a$. 

In some scenarios, we wish to know $\nabla a$ on faces, where $a$
statisfies a Poisson equation; $\nabla^2 a = b$. 

We check if we can design a MG-accelerated iterative solver. Note that
here, the solver *directly* solves for the gradient field.

## Results

~~~gnuplot Convergence
set xr [-0.1: 5.1]
set xlabel 'Cycle nr'
set ylabel 'Max residual'
set grid
set size square
set logscale y
plot 'out' u 1:2 t 'MG Cycle data', 100*0.07**x lw 2 t '0.07^{i}',\
 0.001 lw 2 t 'TOLERANCE'
~~~

Visual inspection

The gradient solver:

![$a_x$ and $a_y$ from the gradient solver](solg/it.png)

Analytical solution:

![$a_x$ and $a_y$ computed from the analytical solution](solg/mg.png)
*/
#include "myQT.h"
#include "utils.h"

scalar b[];
face vector da[];

double l = 2, k = 1;

int main() {
  periodic (left);
  periodic (bottom);
  init_grid (64);
  foreach()
    b[] = -sq(pi)*(sq(2*k) + sq(2*l))*sin(2*k*x*pi)*sin(2*l*y*pi);
  output_ppm (b, file = "b.png", n = 300);

  /**
## The MG-accelerated iterative solver

The solution `da` is iteratively refined with a correction field
`dda`. This is done upto `max_cycles` times.
*/
  double maxres = 0, TOLERANCE = 1e-3; // Residual and its tolerance
  scalar res[], resda[];               // Two residual fields 
  face vector dda[];                   // Solution correction field
  int max_cycles = 100, nrelax = 5;    // Number of MG cycles and relaxation sweeps
  double relaxpar = 1.1;               // Under relaxation paramter
  for (int it = 0; it < max_cycles; it++) {
    /**
       First, the maximum residual is computed on the leaf grid.
    */
    maxres = 0;
    foreach() {
      res[] = b[];
      foreach_dimension() 
	res[] -= (da.x[1] - da.x[])/Delta;
      if (fabs(res[]) > maxres)
	maxres = fabs(res[]);
    }
    printf ("%d %g\n", it, maxres);
    /**
       If the residual does not exceed the tolerance, we can `break`
out of this solver loop.
     */
    if (maxres < TOLERANCE) {
      printf ("# Converged at %d cycles with max_residual = %g\n", it, maxres);
      break;
    }
    /**
### The multigrid V-cycle

The residual is restricted from fine to coarse levels. Then, a loop
starting from the root-cell-level is called. We either take a
zero-initial guess or the face-values are prolongated from the
previous level.
     */
    multigrid_restriction ({res});
    for (int _l = 0; _l <= depth(); _l++) {
      if (_l == 0) {
	foreach_face_level (0) {
	  foreach_dimension()
	    if (is_face_x())
	      dda.x[] = 0;
	}
      } else 
	prolongate_faces_level (dda, _l);
      /**
### Relaxation

The correction field (`dda`) requires its own residual field
(`resda`), using the original `res`-field values. This is needed
because the relaxation itself is independent of neighbor values for
`dda`.
      */
	 for (int rel = 0; rel < nrelax; rel++) {
	   foreach_level(_l) {
	     resda[] = res[];
	     foreach_dimension() 
	       resda[] -= (dda.x[1] - dda.x[])/Delta;
	   }
	   /**
	      A slight under relaxation is advised and can be
	      controlled with `relaxpar`.
	   */
	   foreach_face_level(_l) {
	     foreach_dimension()
	       if (is_face_x())
		 dda.x[] -= Delta*(resda[] - resda[-1])/(2.*dimension*relaxpar);
	   }
	   boundary_level((scalar*){dda}, _l);
	 }
    }
    /**
### Correction
       
The correction of `da` is straighforward.
     */
    foreach_face()
      da.x[] += dda.x[];
    boundary ((scalar*){da});
  } // End of cycles...
  /**
## Visual output
   */
  output_ppm (da.x, file = "daix.png", n = 300);
  output_ppm (da.y, file = "daiy.png", n = 300);
  system ("convert daix.png daiy.png +append it.png");
  vector daa[];
  foreach() {
    daa.x[] = (2*pi*k)*cos(2*pi*k*x)*sin(2*l*pi*y);
    daa.y[] = (2*pi*l)*sin(2*pi*k*x)*cos(2*l*pi*y);
  }
  output_ppm (daa.x, file = "dapx.png", n = 300);
  output_ppm (daa.y, file = "dapy.png", n = 300);
  system ("convert dapx.png dapy.png +append mg.png");
}
