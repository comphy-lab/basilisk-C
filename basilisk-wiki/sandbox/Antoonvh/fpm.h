/**
# Methods for the Finite pointset method

This file implement functions that are relevant for developing methods using the finite pointset method.

We need particles! By default we will assign a  single scalar field "s" to the particles.
 */
#ifndef ADD_PART_MEM
#define ADD_PART_MEM double s;
#endif

#include "particle_reference.h"

/**
   We need a function for finding nearby particles. We will use a helper grid to achieve this efficiently.
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
  foreach_point(X.x, X.y, X.z) {
    foreach_neighbor(neighbor) {
      if (value_p(s, point)) {					
	int * ind = pointer_v(value_p(s, point));			
	for (int n = 0; ind[n] >= 0; n++) {		
	  int _j_particle = ind[n];
	  PARTICLE_VARIABLES;
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
  for (int i = 0; i < nn; i++)
    if (dist[i] < HUGE)
      j++;
  return j;
}
