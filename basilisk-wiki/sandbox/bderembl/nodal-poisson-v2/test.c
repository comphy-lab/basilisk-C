
#define BGHOSTS 2
#include "vertex-utils.h"
#include "nodal-poisson.h"

// qcc -lm -DTRASH=1 -disable-dimensions -grid=multigrid test.c -o test.e

// CC99='mpicc -std=c99' qcc -D_MPI=1  -disable-dimensions -lm test.c -o test.e
// mpirun -np 4 test.e

//GPU on gricad:
//qcc -grid=cuda/multigrid -O2 test.c -o test.e -lm -L$BASILISK/grid/cuda -lbuda -L$BASILISK/grid/gpu -lerrors -lcuda -L/softs/cuda/12.9.patched/lib64 -lnvrtc




vertex scalar omega[];
vertex scalar psi[];

mgstats mgpsi;

int main() {
  N = 1024;
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


  mgpsi = poisson (psi, omega,tolerance=1e-6);
  fprintf(stderr, "i=%d resb=%g resa=%g\n", mgpsi.i, mgpsi.resb, mgpsi.resa);

  // check against exact solution
  double err = 0.;
  foreach_vertex(reduction(max:err)) {
    double psi_exact = -sq(L0)/(5.*sq(pi)) * sin(pi*x/L0) * sin(2.*pi*y/L0);
    err = max(err, fabs(psi[] - psi_exact));
  }
  fprintf(stderr, "max error vs exact: %g\n", err);
}
