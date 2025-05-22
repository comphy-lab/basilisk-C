/**
# Methods for the Finite point(set) method

This file implement functions that are relevant for developing methods
using the finite pointset method.

We need particles! By default we will assign a  single scalar field "s" to the particles.
 */
#ifndef ADD_PART_MEM
#define ADD_PART_MEM double s;
#endif

#include "particle_reference.h"

/**
   We need a function for finding nearby particles. We will use a
   helper grid to achieve this efficiently.
 */
scalar reference[];
bool ref_outdated = true;

// A function that finds the maximum and its index in an array
double max_array (int n, double arr[n], int *i) {
  double max = -HUGE;
  for (int j = 0; j < n; j++) {
    if (arr[j] > max) {
      max = arr[j];
      *i = j;
    }
  }
  return max;
}


int find_nearest_particles (coord X, int nn, Particles plist, int * index, int neighbor = 1) {
  // Update helper field
  if (ref_outdated) {
    free_scalar_data(reference);
    assign_particles (plist, reference);
    ref_outdated = false;
  }
  double dist[nn]; //squared distances;
  for (int i = 0; i < nn; i++)
    dist[i] = HUGE;
  double distm = HUGE; //largest distance
  int il = 0;  //index of largest distance.

  scalar s = reference;
  int _l_particle = plist_s(s);
  // keep a list of nearest particles by forgetting the faw-away ones
  foreach_point(X.x, X.y, X.z) {
    if (point.level > 0) {
      foreach_neighbor(neighbor) {
	foreach_particle_point(reference, point) {
	  double disti = 0;
	  foreach_dimension()
	    disti += sq(X.x - p().x);
	  if (disti < distm) {
	    dist[il] = disti;
	    index[il] = _j_particle;
	    distm = max_array (nn, dist, &il);
	  }
	}
      }
    }
  }
  // count found particles
  int j = 0;
  for (int _i = 0; _i < nn; _i++)
    if (dist[_i] < HUGE)
      j++;
  return j;
}
/**
## Least squares

### p refinement (order refinement)

A lookup table relating the number of particles to the estimation order
*/
int order_barrier_2D[3] = {4, 9, 16}; //2x2, 3x3, 4x4 
int max_particles = 20;
int order_lookup (int n) {
  if (n >= order_barrier_2D[2])
    return 3;
  if (n >= order_barrier_2D[1])
    return 2;
  if (n >= order_barrier_2D[0])
    return 1;
  return 0;
}

// const(1)       s0 = c_0
// linear(3)      s1 = s_0 + c_1x + c_2y  
// quadratic(6)   s2 = s1 + c_3 xy + c_4xx + c_5yy
// cubic(10)      s3 = s2 + c_6 xxy + c_7 xyy + c_8xxx + c_9yyy

int size_lookup (int order) {
  if (order >= 3)
    return 10;
  if (order >= 2)
    return 6;
  if (order >= 1)
    return 3;
  return 1;
}
/**
### The least squares problem matrix

A Vandermonde-like matrix-coefficient finder

We do not apply distance weighing (there is enough to tune already)
 */
#define max_ind (row*m + n)
//Relative to reference location
#define p_x  (pl[part][index[n]].x - loc.x)
#define p_y  (pl[part][index[n]].y - loc.y)
#define p_z  (pl[part][index[n]].z - loc.z)

int fill_matrix_2D (int order = 1, Particles part, int m, int * index, double * A, coord loc) {
  int row = 0;
  if (order >= 0) {
    for (int n = 0; n < m; n++) 
      A[max_ind] = 1.;
  }
  if (order > 0) {
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = p_x;
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = p_y;
  }
  if (order > 1) {
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = p_x*p_y;
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_x);
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_y);
  }
  if  (order > 2) {
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_x)*p_y;
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = p_x*sq(p_y);
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = cube(p_x);
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = cube(p_y);
  }
  return row;
}
/**
### Setup and solve

Using LAPACK
 */
#pragma autolink -lblas -llapack

extern void dgels_(char *trans, int *m, int *n, int *nrhs,
                   double *a, int *lda, double *b, int *ldb,
                   double *work, int *lwork, int *info);

int least_squares_poly_2D (coord loc, double * coefs, Particles parts) {
  int index[max_particles];
  int mat_m =  find_nearest_particles (loc, max_particles, parts, index); // number of equations
  int order = order_lookup (mat_m);
  int mat_n = size_lookup (order); // Number of unknowns
  double * A = malloc(sizeof(double)*mat_m*mat_n);
  fill_matrix_2D (order, parts, mat_m, index, A, loc);
  double rhs[mat_m];
  for (int i = 0; i < mat_m; i++) {
    rhs[i] = pl[parts][index[i]].s;
  }
  // LAPACK stuff suggested by chatGPT, worked on first attempt(!) 
  if (mat_m > 0) { // We a data point
  int lda = mat_m;
  int ldb = mat_m > mat_n ? mat_m : mat_n;
  int info;
  int lwork = -1; // workspace query
  int nrhs = 1;
  double wkopt;
  dgels_("N", &mat_m, &mat_n, &nrhs, A, &lda, rhs, &ldb, &wkopt, &lwork, &info);
  lwork = (int)wkopt;
  double *work = (double *)malloc(lwork * sizeof(double));
  // Actual computation
  dgels_("N", &mat_m, &mat_n, &nrhs, A, &lda, rhs, &ldb, work, &lwork, &info);
  free(work);
  }
  free(A);
  // solution and order 
  for (int i = 0; i < mat_n; i++)
    coefs[i] = rhs[i];
  return order;
}

/**
## smooth operator
 */
void smooth_2D (Particles p) {
  double coefs[9];
  foreach_particle_in(p) {
    coord pc = {x, y};
    int order = least_squares_poly_2D (pc, coefs, p);
    p().z = order > 0 ? coefs[0] : p().s;
  }
  foreach_particle_in(p)
    p().s = p().z;
}

int N_from_part (int np) {
  double snp = sqrt(np);
  int Ni = log(snp)/log(2) - 1;
  return 1<<Ni;
}


/**
## Tests

* [Particle finding](fpm_test.c)
* [Point-based least-sqaures finite difference](fpm_test2.c)
* [Simple advection-equation solver](fpm_advection.c)
* [Covergence test for a diffusion solver](fpm_diffusion.c)

## Todo

* 3D
* Boundaries
* Solve something
* A Poisson solver
* Sensible point distributions (look into packing problem)
* (tree)grid tuning
* combine h and p adaptation
* MPI
* LAPACK in the sandbox
* Look into radial basis functions
*/
