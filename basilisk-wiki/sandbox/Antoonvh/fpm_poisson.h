/**
A solver for the Poisson equation using the finite point method
 */
#ifndef ADD_PART_MEM
#define ADD_PART_MEM double s; double b; double res; double ds; double ds2; 
#endif

#include "fpm.h"
//#include <cblas.h>

int fill_matrix_Poisson_2D (int weight = 1, int order = 1, Particles part, int m, int * index, double * A, coord loc,
		    int * periodic_arr_x = NULL, int * periodic_arr_y = NULL) {
  // No periodicity?
  foreach_dimension() {
    if (periodic_x == false && periodic_arr_x == NULL) 
      periodic_arr_x = (int *)calloc(m, sizeof(int));
  }
  int row = 0;
  int n = m - 1;
  if (order >= 0) {
      for (int n = 0; n < m - 1; n++)  
      A[max_ind] = 1.;
  }
  if (order > 0) {
    row++;
    for (int n = 0; n < m - 1; n++)
      A[max_ind] = p_x;
    row++;
    for (int n = 0; n < m  - 1; n++) 
      A[max_ind] = p_y;
  }
  if (order > 1) {
    row++;
    for (int n = 0; n < m  - 1; n++) 
      A[max_ind] = p_x*p_y;
    row++;
    for (int n = 0; n < m  - 1; n++) 
      A[max_ind] = sq(p_x);
    row++;
    for (int n = 0; n < m  - 1; n++) 
      A[max_ind] = sq(p_y);
  }
  if  (order > 2) {
    row++;
    for (int n = 0; n < m  - 1; n++) 
      A[max_ind] = sq(p_x)*p_y;
    row++;
    for (int n = 0; n < m - 1; n++) 
      A[max_ind] = p_x*sq(p_y);
    row++;
    for (int n = 0; n < m - 1; n++) 
      A[max_ind] = cube(p_x);
    row++;
    for (int n = 0; n < m - 1; n++) 
      A[max_ind] = cube(p_y);
  }
  foreach_dimension()
    if (periodic_x == false && periodic_arr_x == NULL)
      free(periodic_arr_x);
  return row;
}

extern scalar * refl;

double loc_val (coord X, double * coef, int order) {
  
  double val = 0;
  int ind = 0;
  if (order >= 0) {
    val += coef[ind++];
    if (order > 0) {
      val += coef[ind++]*X.x;
      val += coef[ind++]*X.y;
      if (order > 1) {
	val += coef[ind++]*X.x*X.y;
	val += coef[ind++]*sq(X.x);
	val += coef[ind++]*sq(X.y);
	if  (order > 2) {
	  val += coef[ind++]*sq(X.x)*X.y;
	  val += coef[ind++]*X.x*sq(X.y);
	  val += coef[ind++]*cube(X.x);
	  val += coef[ind]*cube(X.y);
	}
      }
    }
  }
  return val;
}

int iterate_particle (Particles parts, int _j_particle, int level = -1) {
  int index[max_particles];
  double dist[max_particles];
  int * periodic_arr_x = NULL;
  int * periodic_arr_y = NULL;
  foreach_dimension() {
    if (periodic_x == true) 
      periodic_arr_x = malloc((max_particles)*sizeof(int));
  }
  scalar s = reference;
  if (level >= 0)
    s = refl[level];
  coord loc = {pl[parts][_j_particle].x, pl[parts][_j_particle].y , 0};
  int neigh = level < 0 ? 1 : 2;
  int mat_m =  find_nearest_particles (loc, max_particles, parts, 
				       index, neigh, periodic_arr_x, periodic_arr_y,
				       true, level = level, reference = s, dist = dist) + 1;
  int ip = 0;
  double md = max_array (mat_m - 1, dist, &ip);
  //fprintf (stderr, "%g %d\n", md, level);
  int order = order_lookup (mat_m - 1);
  int mat_n = size_lookup (order); // Number of unknowns
  double * A = malloc(sizeof(double)*mat_m*mat_n);
  // Weighting of the "Poisson-equation statisfying least squares fit
  double weight = mat_m*sq(md)*1e2;
  int row = fill_matrix_Poisson_2D (weight, order, parts, mat_m, index, A,
				    loc, periodic_arr_x, periodic_arr_y);
  double poisson[10] = {0,0,0,0,2*weight, 2*weight, 0, 0 ,0 ,0};
  int n = mat_m - 1;
  int m = mat_m;
  if (order >= 2) {
    for (int row = 0; row < 6; row++)
      A[max_ind] = poisson[row];
    if (order > 2) {
      for (int row = 6; row < 10; row++)
	A[max_ind] = poisson[row];
    }
  }
  
  double rhs[mat_m];
  rhs[mat_m - 1] = weight*pl[parts][_j_particle].res;
  for (int i = 0; i < mat_m - 1; i++) 
    rhs[i] = pl[parts][index[i]].ds;
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
  /**
Update particles with the Poisson-equation statisfying solution.
   */
  for (int i = 0; i < mat_m -1; i++) {
    coord X;
    foreach_dimension() {
      X.x = pl[parts][index[i]].x - loc.x;
      if (periodic_x)
	X.x += L0*periodic_arr_x[i];
    }
    // Weighing for diagonal domination
    double alphaa = dist[i]/md; 
    pl[parts][index[i]].ds =  alphaa*pl[parts][index[i]].ds + (1 - alphaa)*loc_val(X, rhs, order);
  }
}

double residual (Particles parts) {
  double max_res = 0;
  foreach_particle_in(parts, reduction(max:max_res)) {
    coord X = {x, y, z};
    double coef[10] = {0};
    least_squares_poly_2D(X, coef, parts, true);
    double lap = 2*(coef[4] + coef[5]);
    p().res =  p().b - lap;
    p().ds = 0;
    if (fabs(p().res) > max_res)
      max_res = fabs(p().res);
  }
  foreach_particle_in(parts) {
    double ps = p().s;
    p().s = p().res;
    p().res = ps;
  }
  //smooth_2D (parts);
  foreach_particle_in(parts) {
    double ps = p().s;
    p().s = p().res;
    p().res = ps;
  }
  return max_res;
}

#include "fpm_multi_level.h"

double interpolate_fpm (coord X, Particles parts, int level = -1, scalar s) {
  double coefs[10] = {0};
  int order = least_squares_poly_2D (X, coefs, parts, level = level, reference = s);
  return coefs[0];
}

int sample_particles (Particles p1, Particles p2, int level = -1, scalar s = reference) {
  foreach_particle_in(p1) {
    double temp = p().ds;
    p().ds = p().s;
    p().s = temp;
  }
  foreach_particle_in(p2) {
    coord X = {x,y,z};
    p().ds = interpolate_fpm(X, p1, level, s);
  }
  foreach_particle_in(p1) {
    double temp = p().ds;
    p().ds = p().s;
    p().s = temp;
  }
  
}

int poisson (Particles parts) {
  double resb = residual (parts);
  double TOLER  = 1e-3;
  double resa = resb - 2*TOLER;
  int nrelax = 5;
  int minlevel = 1;
  
  int iter = 0;
  while ((resa > TOLER && (resb - resa) > TOLER) || iter == 0) {
    resb = resa;
    init_multi_level (parts);
    // Multi level cycle
    for (int l = minlevel; l <= depth(); l++) {
      if (l > minlevel) 
	sample_particles (multi_level_parts[l-1], multi_level_parts[l], level = l - 1, s= refl[l-1]);
      for (int i = 0; i < nrelax; i++) {
	foreach_particle_in(multi_level_parts[l]) {
	  iterate_particle(multi_level_parts[l], _j_particle, l);
	}
      }
    }
    // solution particles
    sample_particles (multi_level_parts[depth()], parts, s = refl[depth()]);
    for (int i = 0; i < nrelax; i++) {
      foreach_particle_in(parts) {
	iterate_particle(parts, _j_particle);
      }
    }
    foreach_particle_in(parts) {
      p().s += p().ds;
    }
    resa = residual(parts);
    cleanup_multi_level(parts);
    iter++;
    //printf ("%d %g %g\n", iter, resb, resa);
  }
}
