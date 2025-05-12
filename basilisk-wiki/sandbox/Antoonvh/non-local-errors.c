/**
# Non-local errors for the Poisson equation.

This page concerns the the Poisson solver for $f(x,y)$,

$$\nabla^2 f = s.$$

Since [sine-rich
solutions](http://gerris.dalembert.upmc.fr/gerris/tests/tests/poisson.html) are
*not* suitable for testing, we design another test case where $s$ is
a dipole,

$$s = e^{-(x - 2)^2 - y^2} - e^{-(x + 2.)^2 - y^2}$$ 

The solution ($f$) only exists at the mercy of boundary conditions. 
 */
#include "poisson.h"
#include "utils.h"
scalar s[], f[];

f[top] = dirichlet (0.);
f[bottom] = dirichlet (0.);

#define SOURCE (exp(-sq(x - 2.) - sq(y)) - exp(-sq(x + 2.) - sq(y)))
/**
   The cell-averaged source term is approximated with 64 points.
*/
double source (Point point) {
  double a = 0;
  foreach_child() foreach_child() foreach_child() a += SOURCE;
  return a/pow(pow(2., (double)dimension), 3.);
}
/**
## A reference solution

   Since it may be tricky to find a analytical expression of the
   cell-averaged solution $f$, we compute a reference solution using a
   superior resolution grid. The solution on all relevant levels is
   then stored in an array.
 */

#define OFFSET(l) ((1 << (2*(l)))*((l) + 1))
#define CART_IND(i,j,l) ((i) + (1 << l)*(j))
#define _O (-BGHOSTS - 1)
#define INDEX (OFFSET (level-1) + CART_IND(point.i+_O, point.j+_O, level)) 

int maxlevel = 9;
int reference_level = 11;

double * sol;

void obtain_reference_solution () {
  sol = (double*)malloc (OFFSET (maxlevel)*sizeof(double));
  init_grid (1 << reference_level);
  foreach() 
    s[] = source (point);
  poisson (f, s);
  restriction ({f});
  foreach_cell()
    if (level <= maxlevel)
      sol[INDEX] = f[];
}
/**
## The setup

Periodic conditions are set to exclude the error-inducing effects of
the approximate values of ghost cells. Furthermore we use a large
domain.
*/
int main() {
  periodic (left);
  TOLERANCE = 1e-7;
  L0 = 20.;
  X0 = Y0 = -L0/2;
  
  obtain_reference_solution ();
  /**
We compute a coarse solution.
   */
  init_grid (1 << (maxlevel - 1));
  foreach() 
    s[] = source (point);
  poisson (f, s);
  scalar err_f[];
  foreach()
    err_f[] = f[] - sol[INDEX];
  output_ppm (err_f, n = 512, file = "err.png", min = -1e-5, max = 1e-5);
  /**
     We refine  the center of the domain for a more accurate solution. 
   */
  refine (sq(x) + sq(y) < sq(4) && level < maxlevel);
  foreach() 
    s[] = source (point);
  poisson (f, s);
  foreach() 
    err_f[] = f[] - sol[INDEX];
  output_ppm (err_f, n = 512, file = "err2.png", min = -1e-5, max = 1e-5);
  /**
  Following the suggestion of Steven van der Linden, the oppositely-refined-grid experiment is also performed:
  */
  unrefine (level >= maxlevel - 1);
  refine (sq(x) + sq(y) >= sq(4) && level < maxlevel);
  foreach() 
    s[] = source (point);
  poisson (f, s);
  foreach() 
    err_f[] = f[] - sol[INDEX];
  output_ppm (err_f, n = 512, file = "err3.png", min = -1e-5, max = 1e-5);
}
/**
## Results

The error fields:

![Errors on the relatively coarse equidistant grid](non-local-errors/err.png)

![Errors on the locally refined grid](non-local-errors/err2.png)

![Error on the outer domain refined grid](non-local-errors/err3.png)

Notice that the errors outside the (visibly) refined region decreased
with the refinement elsewhere. This shows that the *source* of the
local error cannot be directly related to a *proper* local error-source
estimate ($\Delta^2 s[0]$).
 */