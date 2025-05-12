/**
# Test for the nodal Poisson solver

![Source](tvpoisson/b.png)

![Solution](tvpoisson/a.png)

Convergence:

~~~bash
Vertex Poisson stats:
mg.i: 5 mg.relax: 5 mg.resa: 7.65397e-07 mg.resb: 0.999715
Cell poisson stats:
mg.i: 6 mg.relax: 3 mg.resa: 6.71873e-06 mg.resb: 0.999805
~~~
 */
#include "nodal-poisson.h"
#include "utils.h"

#define FUNC(a) (exp(-sq(x + (a)) - sq(y - ((a)/2))))
#define SOL (FUNC(-1) - FUNC(2))

int main() {
  L0 = 20;
  X0 = Y0 = -L0/2;
  foreach_dimension()
    periodic (left);
  init_grid (1 << 7);
  refine (fabs(x + 0.5) < 2 && fabs(y - 0.5) < 2 && level <= 8);
  TOLERANCE = 1e-5;
  // Vertex solver:
  scalar a[], b[];
  a.restriction = restriction_vert;
  a.prolongation = refine_vert;
  foreach_vert() {
    a[] = 0;
    b[] = SOL;
  }
  boundary ({a});
  mgstats mg = vpoisson (a, b);
  output_ppm (a, file = "a.png", n = (1 << depth())); 
  output_ppm (b, file = "b.png", n = (1 << depth())); 
  printf("Vertex Poisson stats:\n"
	 "mg.i: %d mg.relax: %d mg.resa: %g mg.resb: %g\n",
	 mg.i, mg.nrelax, mg.resa, mg.resb);
  // Cell solver
  scalar aa[], bb[];
  foreach() {
    aa[] = 0;
    bb[] = SOL;
  }
  boundary ({aa});
  mg = poisson (aa, bb);
  printf("Cell poisson stats:\n"
	 "mg.i: %d mg.relax: %d mg.resa: %g mg.resb: %g\n",
	 mg.i, mg.nrelax, mg.resa, mg.resb);
}
