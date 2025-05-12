/**
# A simple Multigrid Poisson-problem solver

The code in `poisson.h` is very well documented. However, for
modularity, extensability and versatility the code itself is not
structred for the purpose of quicky grasping the miltigrid-accelerated
iterative solver strategy.
*/
#include "grid/multigrid.h"
#include "utils.h"
/**
The discrete Poisson problem can be solved for a single cell by
setting the cell-centered value of a. The idea is that this
"relaxes" the solution a bit towards the correct solution.

~~~literatec
div grad a = b ->
(a[-1] -2*a[] + a[1] + a[0,-1] - 2*a[] + a[0,1] )/sq(Delta) = b[] ->
-4*a[] = b[]*sq(Delta) - (a[-1] + a[1] + a[0,-1] + a[0,1])        ->
a[] = (-b[]*sq(Delta) + (a[-1] + a[1] + a[0,-1] + a[0,1]))/4.;    ->
~~~

The code below relaxes a cell locally. This is called the "direct
replacement method" for relaxation.
 */
void relax_point (Point point, scalar a, scalar b) {
  double ac = -sq(Delta)*b[];
  double d = 2*dimension;
  foreach_dimension()
    ac += (a[-1] + a[1]);
  a[] = ac/d;
}

/** Our solution field `a` does not exactly statisfy the Poisson
equation. We can store the residue in a third scalar field and
assess the maximum absolute redidue.

~~~literatec
residue[] = b[] - (div grad a)[]
~~~

The discrete version of this equation is implemented below.
*/
double residue (scalar a, scalar b, scalar res) {
  double mr = 0; // Maximum residue
  foreach(reduction(max:mr)) {
    res[] = b[];
    foreach_dimension() 
      res[] -= (a[1] + a[-1] - 2*a[])/sq(Delta);
    if (fabs(res[]) > mr)
      mr = fabs(res[]);
  }
  return mr;
}

int main() {
  init_grid (1 << 5);
  /**
We define a dummy problem by setting the field `b` and boundary
conditions for `a`.
   */
  scalar a[], b[];
  periodic(left);
  a[bottom] = dirichlet (0);
  a[top] = dirichlet (0.02);
  foreach() {
    b[] = sin(x/L0*2*pi)*sin(y/L0*2*pi);
    a[] = 0; // Initial guess
  }
  
  /**
The problem will only be solved within a certain tolerance. Further,
in case the iteratively refined solution does not converge, we add an
additional stopping criterion, setting the maximum number of
iterations.
   */
  int MAX_ITERS = 100, it = 0;
  double Tolerance = 1e-4;

  // A field to store the residue
  scalar res[];
  // The main iterative loop, updating the `a` field
  while (residue (a, b, res) > Tolerance && MAX_ITERS > it) {
    /** The residual on coarser levels needs to be computed. This is
	the first leg of the "V" cycle (\), from fine to coarse. Via
	this method, the coarse-level solutions "benifit" from a good
	fine solution, which we will obtain after the first iteration
	and/or from a good initial guess.
    */
    restriction ({res});
    /**
       We will not solve the Poisson problem directly for `a`. Rather
       we will solve for `da` which statisfies the residual equation
       (`div grad da = residue`). Only with this choice can a good
       fine-level solution help with the problem at the coarser
       levels. We also take care that the updated field `da + a`
       statisfies the boundary conditions for the solution field `a`.
     */
    scalar da[];
    for (int b = 0; b < nboundary; b++)
      da.boundary[b] = a.boundary_homogeneous[b];
    /** The second leg of the multigrid cycle "V" cycle (/) iterates from
    coarse fo fine. This way, the fine solutions "benefit" from the
    coarse solutions, as it provides a good initial guess.*/
    for (int l = 0; l <= grid->maxdepth; l++) {
      if (l == 0) // Coarsest guess 
	foreach_level(l)
	  da[] = 0;
      else { // bilinear prolongation of the coarser solution requires boundary conditions
	boundary_level ({da}, l - 1);
	foreach_level(l)
	  da[] = bilinear(point, da);
      }
      // Apply boundary conditions to the guessed solution at this
      // level
      boundary_level ({da}, l);
      /** Relax solution cell-by-cell at this level (single relaxation sweep for simplicity)*/
      foreach_level(l) 
	relax_point (point, da, res);
    } // End of this MG cycle
    /**
The solution field `a` is updated with the `da` field which
approximately statisfies the Poisson equation for the residual.
     */
    foreach() 
      a[] += da[];
    // Log iterative progress for didactic purposes. 
    it++;
    printf ("%d %g\n", it, residue(a, b, res));
  } // end while loop
  
  output_ppm (a, file = "a.png", n = 300, linear = true,
	      map = blue_white_red, min = -0.01, max = 0.03);
  output_ppm (res, file = "res.png", n = 300);
}
/**
## Results

A smooth solution is found

![The solution field](simple_mg/a.png)

It should be noted that the residual (albeit small in magnitude) is
not distributed smoothly.

![The residual field](simple_mg/res.png)

The iterative convergance rate seems quite steady. If could improved
markebly by doing a few more relaxation sweeps on each level.

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
