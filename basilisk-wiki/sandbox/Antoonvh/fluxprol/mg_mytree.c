/**
# A simple Multigrid Poisson-problem solver

*/
#include "myQT.h"
#include "utils.h"

//face vector gar[]; // Gradients or da are not used

void compute_gradients (scalar a, vector ga) {
  foreach_face()
    ga.x[] = face_gradient_x(a, 0);
  flux_prolongate({ga});
}

void relax_point (Point point, scalar a, scalar b) {
  double res = b[];
  double d = 0;
  foreach_dimension() {
#if NOT_SO_CLEVER // Reminder that the effective Neumann boundary
		  // condition on resolution islands is taking the
		  // convergence hostage when it yields an
		  // inconsistent problem at levels.
    if (is_prolongation(neighbor(1)))
      res -= gar.x[1]/(Delta);
    else {
      d++;
      res -= (a[1] - a[0])/sq(Delta);
    } if (is_prolongation(neighbor(-1)))
	res -= -(gar.x[])/(Delta);
    else {
      res -= -(a[] - a[-1])/sq(Delta);
      d++;
    }
#else
    res -= (a[-1] -2*a[] + a[1])/sq(Delta); //uses "dirichlet-like" conditions from prolongation
    d += 2;
#endif
  }
  a[] -= res*sq(Delta)/d;
}

double residue (scalar a, scalar b, scalar res) {
  double mr = 0; // Maximum residue
  face vector gra[];
  compute_gradients(a, gra);
  foreach(reduction(max:mr)) {
    res[] = b[];
    foreach_dimension() 
      res[] -= (gra.x[1] - gra.x[])/Delta;
    if (fabs(res[]) > mr)
      mr = fabs(res[]);
  }
  return mr;
}

int main() {
  init_grid (1 << 5);
  scalar a[], b[];
  a.prolongation = no_restriction;
  periodic(left);
  a[bottom] = dirichlet (0);
  a[top] = dirichlet (0.02);
  foreach() {
    b[] = sin(x/L0*2*pi)*sin(y/L0*2*pi);
    a[] = 0; // Initial guess
  }
  
  int MAX_ITERS = 100, it = 0;
  double Tolerance = 1e-4;

  scalar res[];

  refine (sq(x - 0.4) + sq(y-0.3) < 0.05 && level < 7);
  while (residue (a, b, res) > Tolerance && MAX_ITERS > it) {
    restriction ({res});
    scalar da[];
    for (int b = 0; b < nboundary; b++)
      da.boundary[b] = a.boundary_homogeneous[b];
    for (int l = 0; l <= grid->maxdepth; l++) {
      if (l == 0) // Coarsest guess 
	foreach_level(l)
	  da[] = 0;
      else { // bilinear prolongation of the coarser solution requires boundary conditions
	boundary_level({da}, l - 1);
	foreach_level(l)
	  da[] = bilinear(point, da);
      }
      boundary_level ({da}, l);
      // Sweeping and "level or leaf" does not accelerate convergence
      foreach_level(l)  
	relax_point (point, da, res);
    } // End of this MG cycle
    foreach() 
      a[] += da[];
    it++;
    printf ("%d %g\n", it, residue(a, b, res));
  } // end while loop
  
  output_ppm (a, file = "a.png", n = 300, 
	      map = blue_white_red, min = -0.01, max = 0.03);
  output_ppm (res, file = "res.png", n = 300);
  foreach()
    a[] = level;
  output_ppm (a, file = "lev.png", n = 300);
}
/**
## Results

A smooth solution is found

![The solution field](mg_mytree/a.png)

It should be noted that the residual (albeit small in magnitude) is
not distributed smoothly.

![The residual field](mg_mytree/res.png)

![Refinement level](mg_mytree/lev.png)

The iterative convergence rate seems quite steady. 

~~~gnuplot
set size square
set logscale y
set xr [0:14]
set xlabel 'iteration'
set ylabel 'Maximum residual'
set grid
set key 
plot 'out' t 'data', 10**(-4) t 'Tolerance' lw 2
~~~
*/
