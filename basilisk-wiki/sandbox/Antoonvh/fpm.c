/**
# Finite pointset method

We investigate the idea of a "meshless" approach. Points (particles)
are used as the discrete elements, storing field data alongside their
poisitions.

*/
#define ADD_PART_MEM double s;

#include "particle.h"

# pragma autolink -lblas -llapack
/**
## Linear solver

We need to solve a linear system every time we want to evaluate a
stencil for, e.g., the x-derivative. This routine takes about 1% of
the time-effort in this test.
*/

//LAPACK solver
extern void dgesv_(int * n, int * nrhs, double *a, int *lda,
                   int *ipiv, double *b, int *ldb, int *info);

trace
void solve (int * n, int * nrhs, double *a, int *lda,
	    int *ipiv, double *b, int *ldb, int *info) {
  dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}

/**
## Particle finding

We implement a bad particle-finding algorithm. In this test, this
$\mathcal{O}(n^2)$ routine costs 99% of the total time effort.
 */

trace
void find_nearby_particles (double xp, double yp, Particles fpm, int * index, double suff_d, int needed_points) {
  int i = 0;
  while (i < needed_points) {
    foreach_particle_in(fpm) {
      if (i < needed_points && sq(x -xp) + sq(y - yp) < sq(suff_d) && sq(x -xp) + sq(y - yp) > 0 ) {
	bool to_add = true;
	for (int j = 0; j < i; j++) {
	  if (_j_particle == index[j])
	    to_add = false;
	}
	if (to_add)
	  index[i++] = _j_particle;
      }
    }
    suff_d *= 1.2;
  }
}

int main() {
  // A helper grid which is not used
  L0 = 5;
  X0 = Y0 = -L0/2.;
  init_grid (N);

  FILE * fp = fopen ("data", "w");
  // Convergence test
  for (int npart = 100; npart < 2e5; npart *= 2) {
    // Initialize particles with gaussian blob data
    Particles fpm = new_particles (npart);
    foreach_particle_in(fpm) {
      foreach_dimension()
	p().x = L0*noise()/2;
    }
    foreach_particle_in(fpm) {
      p().s = exp(-sq(x) - sq(y));
    }
    /**
       ~~~gnuplot Initialized pointset
       set size square
       plot 'data' u 1:2:3 palete pt 5
       ~~~
       
       We will now try to estimate the x-derivative at each point using neighboring points.
    */
    int needed_points = 3; // for 3 unknowns s = ax + by + c 
    double suff_distance = L0/sqrt(npart);
    
    double error = 0;
    foreach_particle_in(fpm, reduction(+:error)) {
      int index[needed_points];
      find_nearby_particles (x, y, fpm, index, suff_distance, needed_points - 1);
      int n = 3;        // Number of equations
      int nrhs = 1;     // Number of right-hand sides
      int lda = 3;      // Leading dimension of a
      int ldb = 3;      // Leading dimension of b
      int info;
      int ipiv[3];      // Pivot locations
      // rhs vector
      double b[3] = {p().s,  pl[fpm][index[0]].s,  pl[fpm][index[1]].s};
      // matrix appears transposed
      double A[9] = {0, pl[fpm][index[0]].x - x, pl[fpm][index[1]].x - x,
		     0, pl[fpm][index[0]].y - y, pl[fpm][index[1]].y - y,
		     1, 1, 1};
      solve(&n, &nrhs, A, &lda, ipiv, b, &ldb, &info);
      if (info != 0)
	fprintf (stderr, "LAPACK Solver error \n");

      if (npart == 3200) 
	fprintf (fp, "%g %g %g %g %g %g\n", x, y, p().s, b[0], b[1], b[2]);
      
      error += fabs(b[0] + 2*x*exp(-sq(x) - sq(y)));
    }
    printf ("%d %g\n",npart, error/npart);
  }
  fclose (fp);
}

/**
## Results

 ~~~gnuplot Estimated x-derivative 
       set size square
       set sbrange [-1:1]
       plot 'data' u 1:2:4 palete pt 5
 ~~~

 ~~~gnuplot Error convergence is first order? 
 set size square
 set logscale xy
 set grid
 set xlabel 'sqrt(N)'
 set ylabel 'Mean error'
 plot 'out' u ($1**(0.5)):2, 10*x**(-1)
 ~~~
*/
 
      
