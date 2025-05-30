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
//periodic switches
bool periodic_x = false;
bool periodic_y = false;
bool periodic_z = false;
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

/**
## Finding nearby particles

We will lookup nearby particles that are positioned in neighboring
cells in the helper grid. In order to increase our search distance, we
can may look at coarser levels.
 */
Point locate_level (double xp = 0., double yp = 0., double zp = 0., int level = -1){
  for (int l = (level >= 0 ? level : depth()); l >= 0; l--) {
    Point point = {0};
    point.level = l;
    int n = 1 << point.level;
    point.i = (xp - X0)/L0*n + GHOSTS;
#if dimension >= 2
    point.j = (yp - Y0)/L0*n + GHOSTS;
#endif
#if dimension >= 3
    point.k = (zp - Z0)/L0*n + GHOSTS;
#endif
    if (point.i >= 0 && point.i < n + 2*GHOSTS
#if dimension >= 2
	&& point.j >= 0 && point.j < n + 2*GHOSTS
#endif
#if dimension >= 3
	&& point.k >= 0 && point.k < n + 2*GHOSTS
#endif
	) {
      if (allocated(0) && is_local(cell) && (is_leaf(cell) || is_coarse()))
      return point;
    }
   else break;
  }
  Point point = {0};
  point.level = -1;
  return point;
}

macro2 foreach_point_level (double _x = 0., double _y = 0., double _z = 0.,
			    int level = -1, char flags = 0, Reduce reductions = None)
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    coord _p = { _x, _y, _z };
    Point point = locate_level (_p.x, _p.y, _p.z, level); // fixme
    if (point.level >= 0)
      {...}
  }
}

/**
The function that finds (the indices) of nearby points
 */

int find_nearest_particles (coord X, int nn, Particles plist, int * index, 
			    int neighbor = 1, int * index_periodic_x = NULL,
			    int * index_periodic_y = NULL, bool self = true,
			    int level = -1, scalar reference = reference,
			    double * dist = NULL) {
  // Take care of periodically placed particles
  coord dim = {0,1,2};
  foreach_dimension() {
    if (index_periodic_x == NULL && periodic_x == true)
      fprintf (stderr, "Periodic in dim %d, bit no alloctated indices array\n", dim.x);
    if (periodic_x == true) {
      for (int i = 0; i < nn; i++)
	index_periodic_x[i] = 0;
    }
  }
  // Update helper field
  if (ref_outdated) {
    free_scalar_data(reference);
    assign_particles (plist, reference);
    ref_outdated = false;
  }
  bool alloc = false;
  if (dist == NULL) {
    alloc = true;
    dist = malloc(nn*sizeof(double)); //squared distances;
  }
  for (int i = 0; i < nn; i++)
    dist[i] = HUGE;
  double distm = HUGE; //largest distance
  int il = 0;  //index of largest distance.
  
  scalar s = reference;
  int _l_particle = plist_s(s);
  // keep a list of nearest particles by forgetting the faw-away ones, i.e. no sorting.
  foreach_point_level(X.x, X.y, X.z, level) {
    if (point.level > 0) {
      //printf ("new point\n");
      foreach_neighbor(neighbor) {
	//printf ("_k = %d, _l = %d\n", _k, _l);
	coord ni = {_k, _l, 0};
	//check_double (point, reference);
	foreach_particle_point(reference, point) {
	  double disti = 0;
	  foreach_dimension() {
	    double px = p().x;
	    if (periodic_x) {
	      if (ni.x > 0 && p().x < X.x)
		px += L0;
	      else if (ni.x < 0 && p().x > X.x)
		px -= L0;
	    }
	    disti += sq(X.x - px);
	  }
	  if (disti < distm && disti < 2*sq(neighbor*Delta) && (self || disti > Delta*1e-9)) {
	    dist[il] = disti;
	    index[il] = _j_particle;
	    foreach_dimension() {
	      if (periodic_x) {
		index_periodic_x[il] = 0;
		if (ni.x > 0 && p().x < X.x)
		  index_periodic_x[il] = 1;
		else if (ni.x < 0 && p().x > X.x)
		  index_periodic_x[il] = -1;
	      }
	    }
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
  if (alloc)
    free(dist);
  return j;
}
/**
## Least squares

### p refinement (order refinement)

A lookup table relating the number of particles to the estimation order (-1)
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

#define p_x  (pl[part][index[n]].x - loc.x + L0*periodic_arr_x[n])
#define p_y  (pl[part][index[n]].y - loc.y + L0*periodic_arr_y[n])
#define p_z  (pl[part][index[n]].z - loc.z + L0*periodic_arr_z[n])

int fill_matrix_2D (int order = 1, Particles part, int m, int * index, double * A, coord loc,
		    int * periodic_arr_x = NULL, int * periodic_arr_y = NULL) {
  // No periodicity?
  foreach_dimension() {
    if (periodic_x == false && periodic_arr_x == NULL) 
      periodic_arr_x = (int *)calloc(m, sizeof(int));
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
  foreach_dimension()
    if (periodic_x == false && periodic_arr_x == NULL)
      free(periodic_arr_x);
  return row;
}
/**
### Setup and solve

Using LAPACK
 */
#pragma autolink -lcblas -llapack

extern void dgels_(char *trans, int *m, int *n, int *nrhs,
                   double *a, int *lda, double *b, int *ldb,
                   double *work, int *lwork, int *info);

int least_squares_poly_2D (coord loc, double * coefs, Particles parts,
			   bool self = true, int level = -1, scalar reference = reference) {
  int index[max_particles];
  int * periodic_arr_x = NULL;
  int * periodic_arr_y = NULL;

  foreach_dimension() {
    if (periodic_x == true) 
      periodic_arr_x = malloc((max_particles)*sizeof(int));
  }
   // number of equations mat_m
  int mat_m =  find_nearest_particles (loc, max_particles, parts, 
				       index, 1, periodic_arr_x, periodic_arr_y,
				       self, level, reference);
  int order = order_lookup (mat_m);
  int mat_n = size_lookup (order); // Number of unknowns
  double * A = malloc(sizeof(double)*mat_m*mat_n);
      
  fill_matrix_2D (order, parts, mat_m, index, A, loc, periodic_arr_x, periodic_arr_y);
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
  for (int i = 0; i < mat_n; i++)
    coefs[i] = rhs[i];

  foreach_dimension() {
    if (periodic_x == true) 
      free(periodic_arr_x);
  }
  
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
