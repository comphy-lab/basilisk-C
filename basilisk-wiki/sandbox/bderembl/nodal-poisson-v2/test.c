#include "grid/multigrid.h"
#include "vertex-utils.h"
#include "nodal-poisson.h"

// qcc -lm -DTRASH=1 -disable-dimensions test.c -o test.e

// CC99='mpicc -std=c99' qcc -D_MPI=1  -disable-dimensions -lm test.c -o test.e
// mpirun -np 2 test.e






vertex scalar omega[];
vertex scalar psi[];

mgstats mgpsi;

int main() {
  N = 8;
  // if run in parallel
//  dimensions(nx=2, ny=1);

  init_grid (N);

  omega[left]   = 0.;
  omega[right]  = 0.;
  omega[top]    = 0.;
  omega[bottom] = 0.;

  psi[left]   = 0.;
  psi[right]  = 0.;
  psi[top]    = 0.;
  psi[bottom] = 0.;


  double omega0 = 1.0 [0,-1];
  double psi0   = 0.0 [2,-1];

  foreach_vertex() {
    psi[]   = psi0;
    omega[] = omega0*sin(pi*x/L0)*sin(2*pi*y/L0);
  }

  boundary({psi, omega});

  /* foreach_vertex() { */
  /*   printf("%g \t %g \t %g \t %g\n",x,y,psi[],omega[]); */
  /* } */


  mgpsi = poisson (psi, omega);

}
