/**

# How to solve $\psi_{xx} - c\psi = 0$, with Baslisk. 

The aforementioned equation is a version of the 1D Poisson-Helmholtz
(PH) equation:

$$L(a) = \nabla\cdot (\alpha\nabla a) + \lambda a = b,$$

and hence we include: 
*/

#include "grid/multigrid1D.h"
#include "poisson.h"

/**
Scalar fields to store our solution and the
right-hand-side of the full PH equation ($b$) are initialized.
*/
scalar psi[], s[]; 
/**
And boundary conditions are set to arive at a non-trivial solution.
*/
psi[left] = dirichlet(1);
psi[right] = neumann(0.);
/**
The equation is solved in the `main()` function, using 128 cells;
*/
int main() {
  init_grid (128);
  /**
  We choose $c = 2$ and solve the equation via the `poisson()` user
  interface function.
  */
  const scalar c[] = -2.;
  poisson (psi, s, lambda = c);
  /**
Finally, the solution is printed so that it may be plotted:
  */ 
  foreach()
    printf("%g %g\n", x, psi[]);
}
/**
Resulting in, 

~~~gnuplot the result.
set xlabel 'x'
set ylabel '{/Symbol psi}'
set key off
set grid
set size ratio -1
plot 'out' w l lw 3
~~~
 */
