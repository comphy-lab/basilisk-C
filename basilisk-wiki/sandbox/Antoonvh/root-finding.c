/**
# Root Finding Using Adaptive Grids

We can find the zeros of a smooth function using adaptive grids. This
can be achieved by iteratively refining arround the neighborhood where
the function changes sign. 

In this example we try to find the non trivial zeros of the first
Bessel function of the first kind for $x$ is between zero and 20. I.e.

$$J_1(x) = 0,\ \mathrm{for  }\ x \in(0,20).$$ 

We assume that the `math` library has a sufficiently accurate definitnion of
this function.
*/

#include "grid/bitree.h"

#define FUN(x) (j1(x))
scalar near[]; 
int maxlevel = 19, baselevel = 3; 

int main() {
  L0 = 20;
  X0 = 1e-9;
  init_grid (1 << baselevel);
  // Iteratively refine the grid
  for (int l = baselevel; l < maxlevel; l++) {
    foreach_level (l) {
      near[] = 0; 
      if (FUN(x - Delta/2.)*FUN(x + Delta/2.) <= 0)
	near[] = (double)true;
    }
    refine (near[] == (double)true && level <= l);
  }
  // Output the zeros:
  foreach_level (depth())
    if (FUN(x - Delta/2.)*FUN(x + Delta/2.) <= 0)
      printf ("%.8g\n", x);
  printf ("Error margin: +/- %.2g\n", L0/(1 << depth()));
  // Count the used cells (including parents)
  int n = 0;
  foreach_cell()
    n++;
  printf ("Number of Cells: %d\n"
	  "Equdistant-grid equivalent number: %d\n", n, 1 << depth());
}

/**
## Results

This produces the following output in the terminal:

~~~bash
3.8316917
7.0155907
10.173473
13.323689
16.470623
19.615841
Error margin: +/- 3.8e-05
Number of Cells: 853
Equdistant-grid equivalent number: 524288
~~~

This seems to correspond with online references down to the mentioned accuracy.  
*/
