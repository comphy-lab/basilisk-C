/**
## Reference particles interface for fpm


 */

/**
We need particles! By default we will assign a  single scalar field "s" to the particles.
 */
#ifndef ADD_PART_MEM
#define ADD_PART_MEM double s; double * sl;
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
			    int * index_periodic_y = NULL, int * index_periodic_z = NULL, bool self = true,
			    int level = -1, scalar reference = reference,
			    double * dist = NULL) {
  // Take care of periodically placed particles
  coord dim = {0,1,2};
  foreach_dimension() {
    if (index_periodic_x == NULL && periodic_x == true)
      fprintf (stderr, "Periodic in dim %g, bit no alloctated indices array\n", dim.x);
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
  NOT_UNUSED(_l_particle);
  // keep a list of nearest particles by forgetting the faw-away ones, i.e. no sorting.
  foreach_point_level(X.x, X.y, X.z, level) {
    if (point.level > 0) {
      //printf ("new point\n");
      foreach_neighbor(neighbor) {
	//printf ("_k = %d, _l = %d\n", _k, _l);
#if dimension == 2
	coord ni = {_k, _l, 0};
#else // 3 dims
	coord ni = {_l, _m, _n};
#endif

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
	  if (disti < distm && disti < dimension*sq((neighbor + 0.5)*Delta) && (self || disti > Delta*1e-9)) {
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
