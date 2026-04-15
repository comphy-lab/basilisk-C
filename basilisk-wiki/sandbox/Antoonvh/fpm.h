/**
# Methods for the Finite point(set) method

This file implement functions that are relevant for developing methods
using the finite pointset method.

We need some functions to localize particles
*/
#include "fpm_ref.h"
/**
## Least squares

### p refinement (order refinement)

A lookup table relating the number of particles to the estimation order (-1)
*/
# if dimension == 2
int order_barrier[3] = {4, 9, 16}; //2x2, 3x3, 4x4
# else 
int order_barrier[3] = {8, 18, 30};
#endif
int max_particles = 20;
int order_lookup (int n) {
  if (n >= order_barrier[2])
    return 3;
  if (n >= order_barrier[1])
    return 2;
  if (n >= order_barrier[0])
    return 1;
  return 0;
}

// const(1)       s0 = c_0
// linear(3)      s1 = s_0 + c_1x + c_2y  
// quadratic(6)   s2 = s1 + c_3 xy + c_4xx + c_5yy
// cubic(10)      s3 = s2 + c_6 xxy + c_7 xyy + c_8xxx + c_9yyy
// 3D
// linear(4)     s1 = c_0 + c_1x + c_2y + c_3z  
// quadratic(10)   s2 = s1 + c_4 xy + c_5 xz + c_6 yz + c_7 xx + c_8 yy + c_9 zz
// cubic(20)      s3 = s2 + c_10 xxy + c_11 xyy + c_12 xxz + c_13 xzz + c_14 yzz + c_15 yyz + c_16 xyz + c_17xxx + c_18yyy + c_19 zzz


int size_lookup (int order) {
  if (order >= 3)
    return dimension == 2 ? 10 : 20;
  if (order >= 2)
    return dimension == 2 ? 6 : 10;
  if (order >= 1)
    return dimension == 2 ? 3 : 4;
  return 1;
}
/**
### The least squares problem matrix

A Vandermonde-like matrix-coefficient finder

We do not apply distance weighing (there is enough to tune already)
 */
#define max_ind (row*m + n)
//Relative to reference location

#define p_x  (pl[part][index[n]].x - loc.x + L0*periodic_arr_x[n])
#define p_y  (pl[part][index[n]].y - loc.y + L0*periodic_arr_y[n])
#define p_z  (pl[part][index[n]].z - loc.z + L0*periodic_arr_z[n])

int fill_matrix (int order = 1, Particles part, int m, int * index, double * A, coord loc,
		    int * periodic_arr_x = NULL, int * periodic_arr_y = NULL, int * periodic_arr_z = NULL) {
  // No periodicity?
  bool alloc_x = false;
  bool alloc_y = false;
  bool alloc_z = false;
  foreach_dimension() {
    if (periodic_x == false && periodic_arr_x == NULL) {
      alloc_x = true;
      periodic_arr_x = (int *)calloc(m, sizeof(int));
    }
  }
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
#if dimension > 2
    row++;
    for (int n = 0; n < m; n++)
      A[max_ind] = p_z;
#endif
  }
  if (order > 1) {
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = p_x*p_y;
#if dimension > 2
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = p_x*p_z;
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = p_y*p_z;
#endif
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_x);
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_y);
#if dimension > 2
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_z);
#endif
  }
  if  (order > 2) {
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_x)*p_y;
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = p_x*sq(p_y);
#if dimension > 2
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_x)*p_z;
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_z)*p_x;
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_y)*p_z;
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = sq(p_z)*p_y;
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = p_x*p_y*p_z;
#endif
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = cube(p_x);
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = cube(p_y);
#if dimension > 2
    row++;
    for (int n = 0; n < m; n++) 
      A[max_ind] = cube(p_z);
#endif
  }
  foreach_dimension()
    if (alloc_x == true)
      free(periodic_arr_x);
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

int least_squares_poly (coord loc, double * coefs, Particles parts,
			bool self = true, int level = -1, scalar reference = reference,
			int * list = NULL) {
  int index[max_particles];
  int * periodic_arr_x = NULL;
  int * periodic_arr_y = NULL;
  int * periodic_arr_z = NULL;
  foreach_dimension() {
    if (periodic_x == true) 
      periodic_arr_x = malloc((max_particles)*sizeof(int));
  }
   // number of equations mat_m
  int mat_m =  find_nearest_particles (loc, max_particles, parts, 
				       index, 1, periodic_arr_x, periodic_arr_y,
				       periodic_arr_z, self, level, reference);
  //printf ("%d\n", mat_m);
  int order = order_lookup (mat_m);
  int mat_n = size_lookup (order); // Number of unknowns
  double * A = malloc(sizeof(double)*mat_m*mat_n);
      
  fill_matrix (order, parts, mat_m, index, A, loc, periodic_arr_x, periodic_arr_y, periodic_arr_z);
  int c = 1;
  if (list != NULL) {
    c = 0;
    while (list[c++] >= 0);
    c--;
  }
  double rhs[c*mat_m];
  for (int j = 0; j < c; j++) 
    for (int i = 0; i < mat_m; i++) {
      if (list == NULL) {
	rhs[i] = pl[parts][index[i]].s;
	//printf ("rhs = %g\n", rhs[i]);
      }
	else
	rhs[j*mat_m + i] = pl[parts][index[i]].sl[list[j]];
    }
  
  // LAPACK stuff suggested by chatGPT, worked on first attempt(!) 
  int ldb = mat_m > mat_n ? mat_m : mat_n;
  //printf ("%d %d %d\n", mat_m, mat_n, c);
  if (mat_m > 0) { // We a data point
    int lda = mat_m;
    int info;
    int lwork = -1; // workspace query
    int nrhs = c;
    double wkopt;
    // printf ("%d %d\n", mat_m, mat_n);
    
    dgels_("N", &mat_m, &mat_n, &nrhs, A, &lda, rhs, &ldb, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    double *work = (double *)malloc(lwork * sizeof(double));
    // Actual computation
    // printf ("a %d %d\n", mat_m, mat_n);
    dgels_("N", &mat_m, &mat_n, &nrhs, A, &lda, rhs, &ldb, work, &lwork, &info);
    free(work);
  }
  free(A);
  // solution and order 
  for (int j = 0; j < c; j++)
    for (int i = 0; i < ldb; i++) {
      coefs[j*mat_n + i] = rhs[i + j*ldb];
    }
  foreach_dimension() {
    if (periodic_x == true) 
      free(periodic_arr_x);
  }
  return mat_n;
}

/**
## smooth operator
 */
void smooth_2D (Particles p) {
  double coefs[9];
  foreach_particle_in(p) {
    coord pc = {x, y};
    int order = least_squares_poly (pc, coefs, p);
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

event defaults (i = 0) {
  foreach_dimension() {
    reference[left] = 0;
    reference[right] = 0;
  }
  reference.restriction = ref_restriction;
  reference.prolongation = ref_prolongation;
#if TREE
  reference.coarsen = ref_coarsen;
  reference.refine = ref_prolongation;
#endif
}


void free_p_fpm (Particles p) {
  foreach_particle() {
    if (p().sl != NULL)
      free(p().sl);
  }
  free_p();
}

event cleanup (t = end) {
    free_scalar_data (reference);
}



/**
## Tests

* [Particle finding](fpm_test.c)
* [Point-based least-sqaures finite difference](fpm_test2.c)
* [Simple periodic advection-equation solver](fpm_advection.c)
* [Covergence test for a diffusion-equation solver](fpm_diffusion.c)
* [Multi level particle finding](fpm_multi_level.c)

## Todo

* 3D
* Boundaries
  - ~Periodic box boundaries~
  - Arbitrary Dirichlet boundary particles
  - Arbitrary Von Neumann boundary particles
* Solve something
* A Poisson solver
* A Lagrangian 2D vorticity-equation solver
* Sensible point distributions (look into packing problem)
* (tree)grid tuning
* combine h and p adaptation
* MPI
* fopenmp
* LAPACK in the sandbox
* Look into radial basis functions
*/
