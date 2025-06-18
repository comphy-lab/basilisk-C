/**
A solver for the Poisson equation using the finite point method
 */
#ifndef ADD_PART_MEM
#define ADD_PART_MEM double s; double b; double res; double ds; 
#endif

#include "fpm.h"
#include "fpm_double_ls.h"

extern scalar * refl;

int n_smallest_of_m (int m, double a[m], int index[m], int n, int indexo[n]) {
  double ao[n];
  for (int i = 0; i < n; i++)
    ao[i] = HUGE;
  double largest = HUGE;
  int ind = 0;
  for (int i = 0; i < m; i++) {
    if (a[i] < largest) {
      ao[ind] = a[i];
      indexo[ind] = index[i];
      largest = max_array (n, ao, &ind);
    }
  }
  int j = 0;
  for (int i = 0; i < n; i++)
    if (ao[i] < HUGE)
      j++;
  return j;
}


//P is an np*10 matrix, in columnform, computing the estimated laplacian.
void Laplacian_matrix (coord loc, Particles part,  int np, int * index, double * P, int * periodic_arr_x, int * periodic_arr_y, double * dis, double md) {
  
  bool alloc_x = false;
  bool alloc_y = false;
  
  foreach_dimension() {
    if (periodic_x == false) {
      alloc_x = true;
      periodic_arr_x = (int *)calloc(np, sizeof(int));
    }
  }
  
  for (int n = 0; n < np; n++) {
    int j = 0;
    P[n + np*j++] = 0; //c_0 a -> 0
    P[n + np*j++] = 0; //c_1 x -> 0
    P[n + np*j++] = 0; //c_2 y -> 0
    P[n + np*j++] = 0; //c_3 xy -> 0
    P[n + np*j++] = (1 - dis[n]/md)*2.; //c_4 xx -> 2
    P[n + np*j++] = (1 - dis[n]/md)*2.; //c_5 yy -> 2
    P[n + np*j++] = (1 - dis[n]/md)*2.*p_y; //c_6 xxy -> 2y
    P[n + np*j++] = (1 - dis[n]/md)*2.*p_x; //c_7 yyx -> 2x
    P[n + np*j++] = (1 - dis[n]/md)*6.*p_x; //c_8 xxx -> 6x
    P[n + np*j]   = (1 - dis[n]/md)*6.*p_y; //c_9 yyy -> 6y
  }
  
  foreach_dimension() {
    if (alloc_x) {
      free(periodic_arr_x);
    }
  }
}



int iterate_particle (Particles parts, int _j_particle, int level = -1) {
  int index[max_particles];
  double dist[max_particles];
  int * periodic_arr_x = NULL;
  int * periodic_arr_y = NULL;
  int d = 0;
  foreach_dimension() {
    if (periodic_x == true) {
      periodic_arr_x = malloc((max_particles)*sizeof(int));
    }
  }
  scalar s = reference;
  if (level >= 0)
    s = refl[level];
  coord loc = {pl[parts][_j_particle].x, pl[parts][_j_particle].y , 0};
  int neigh = level < 0 ? 1 : 2;
  int mat_m =  find_nearest_particles (loc, max_particles, parts, 
				       index, neigh, periodic_arr_x, periodic_arr_y,
				       true, level = level, reference = s, dist = dist) ;
  
  int ip = 0;
  double md = max_array (mat_m , dist, &ip);
  //fprintf (stderr, "%g %d\n", md, level);
  int order = order_lookup (mat_m );
  int mat_n = size_lookup (order); // Number of unknowns
  double * A = malloc(sizeof(double)*mat_m*mat_n);
  int row = fill_matrix_2D (order, parts, mat_m, index, A,
			    loc, periodic_arr_x, periodic_arr_y);
  int nc = 1;//mat_m;
  int indc[nc];
  nc =  n_smallest_of_m (mat_m, dist, index, nc, indc);
  int period_x[nc];
  int period_y[nc];
  double dis[nc];
  for (int i = 0; i < nc; i++) {
    for (int j = 0 ; j < mat_m; j++)
      if (indc[i] == index[j]) {
	period_x[i] = periodic_arr_x[j];
	period_y[i] = periodic_arr_y[j];
	dis[i] = dist[j];
      }
  }
  // printf ("%d %d %d %d %d\n", _j_particle, indc[0], order, nc, mat_m);
  double poisson[nc*10];// = malloc(nc*10*sizeof(double));
  double a[nc];
  Laplacian_matrix (loc, parts, nc , indc, poisson, period_x, period_y, dis, md);
  for (int i = 0; i < nc; i++)
    a[i] = pl[parts][indc[i]].res*(1 - dis[i]/md);
  int n = mat_n;
  int m = mat_m;
  double ds[m];
  for (int j = 0; j < m; j++)
    ds[j] = pl[parts][index[j]].ds;
  update_double_ls (A, m, n, ds, poisson, nc, a);
  free (A);
  /**
Update particles with the Poisson-equation statisfying solution.
   */
  for (int i = 0; i < mat_m; i++) {
#if FPM_BOUNDARY
    // Skipping boundary points in the neighborhood
    if (pl[parts][index[i]].bound == 0) {
#endif
      double alphaa = 1;//(1 - dist[i]/md)/mat_m; 
      pl[parts][index[i]].ds += alphaa* ds[i];
#if FPM_BOUNDARY
    }
#endif
  }
  foreach_dimension() {
    if (periodic_x == true && periodic_arr_x != NULL) {
      free(periodic_arr_x);
    }
  }
   
}

double residual (Particles parts) {
  double max_res = 0;
  foreach_particle_in (parts, reduction(max:max_res)) {
    p().ds = 0;
#if FPM_BOUNDARY
    // Skipping boundary points 
    if (p().bound == 0) {
#endif
      coord X = {x, y, z};
      double coef[10] = {0};
      least_squares_poly_2D(X, coef, parts, true);
      double lap = 2*(coef[4] + coef[5]);
      p().res =  p().b - lap;
      if (fabs(p().res) > max_res)
	max_res = fabs(p().res);
#if FPM_BOUNDARY
    }
#endif
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
#if FPM_BOUNDARY
    // Skipping boundary points 
    if (p().bound == 0) {
#endif
    coord X = {x,y,z};
    p().ds = interpolate_fpm(X, p1, level, s);
#if FPM_BOUNDARY
    }
#endif
  }
  foreach_particle_in(p1) {
    double temp = p().ds;
    p().ds = p().s;
    p().s = temp;
  }
}

double TOLER  = 1e-3;

int poisson (Particles parts) {
  double resb = residual (parts);
  printf ("resb = %g\n", resb);
  double resa = resb - 2*TOLER;
  int nrelax = 4;
  int minlevel = 2;
  
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
#if FPM_BOUNDARY
	  // Skipping boundary points 
	  if (p().bound == 0) {
#endif
	  iterate_particle(multi_level_parts[l], _j_particle, l);
#if FPM_BOUNDARY
	  }
#endif
	}
      }
    }
    // solution particles
    sample_particles (multi_level_parts[depth()], parts, s = refl[depth()]);
    for (int i = 0; i < nrelax; i++) {
      foreach_particle_in(parts) {
#if FPM_BOUNDARY
	// Skipping boundary points 
	if (p().bound == 0) {
#endif
	iterate_particle(parts, _j_particle);
#if FPM_BOUNDARY
	  }
#endif
      }
    }
    double max_ds = 0;
    foreach_particle_in(parts) {
      p().s += p().ds;
      if (fabs(p().ds) > max_ds)
	max_ds = fabs(p().ds);
    }
    resa = residual(parts);
    cleanup_multi_level(parts);
    iter++;
    printf ("%d %g %g %g\n", iter, resb, resa, max_ds);
  }
}


