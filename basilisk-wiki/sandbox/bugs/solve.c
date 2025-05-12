/**
# `SIGSEGV` in  `solve()`

A minimal example.  

In other versions , `nan`s arise. 

*/
#include "grid/multigrid.h"
#include "solve.h"

void solver (scalar b, vector v) {
  scalar not_used[];
  scalar da = v.x;
  solve (da, (da[0,-1] + da[0,1])/10. + da[], b[]);
}

int main() {
  N = 16;
  init_grid (N);
  scalar non_used2[], b[];
  vector v[];
  foreach() {
    b[] = 0;
    foreach_dimension()
      v.x[] = 0;
  }
  solver(b, v);
}
