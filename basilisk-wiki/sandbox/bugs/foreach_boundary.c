#include "grid/tree.h"
#include "run.h"

scalar T[];

int main() {
  L0 = 1;
  init_grid (128);
  run();
}

T[top] = dirichlet(-0.5);

event init(t=0){
  foreach()
    T[] = 1.;
  boundary({T});
}

event flux(t=0){
  foreach_boundary (top)
    fprintf(stdout, "%g %g", T[ghost], T[]);
}