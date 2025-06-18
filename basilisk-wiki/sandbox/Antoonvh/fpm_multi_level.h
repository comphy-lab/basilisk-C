/**
# Finite point Multi level stuff

To aid with the linear solver.
 */

int * multi_level_parts;

#ifndef ADD_PART_MEM
#define ADD_PART_MEM double s; double b; double res; double ds; double ds2; 
#endif
#include "particle_reference.h"
scalar * refl = NULL;

int collect_particles_point (Point point, Particles p, Particles new_list) {
  int np = 0;
  particle new_p = {.x = 0, .y = 0, .z = 0, .s = 0, .res = 0, .ds = 0}; 
  foreach_particle_point(reference, point) {
#if FPM_BOUNDARY
    if (p().bound == 0) {
#endif
      np++;
      foreach_dimension()
	new_p.x += p().x;
      new_p.res += p().res;
#if FPM_BOUNDARY
    }
#endif
  }
  if (np) {
    foreach_dimension()
      new_p.x /= np;
    new_p.res /= np;
    new_p.ds  = 0;
    add_particle(new_p, new_list);
  }
#if FPM_BOUNDARY
  // dirichlet particles
  {
    np = 0;
    new_p = (particle){.x = 0, .y = 0, .z = 0, .s = 0, .res = 0, .ds = 0, .bound = 1};
    foreach_particle_point(reference, point) {
      if (p().bound == 1) {
	np++;
	foreach_dimension()
	  new_p.x += p().x;
	new_p.res += p().res;
	new_p.s += p().s;
      }
    }
  if (np) {
    foreach_dimension()
      new_p.x /= np;
    new_p.res /= np;
    new_p.ds  = 0;
    new_p.s /= np;
    add_particle(new_p, new_list);
  }
  }
#endif
  return np;
}

int init_multi_level(Particles p) {
  if (ref_outdated) {
    assign_particles (p, reference);
    ref_outdated = false;
  }
  int dep = grid->maxdepth;
  multi_level_parts = malloc (sizeof(int)*(dep + 1));
  for (int l = 0; l <= dep; l++) {
    scalar s = new_scalar("ref");
    refl = list_add (refl, s);
    multi_level_parts[l] = new_particles(0);
    foreach_level_or_leaf(l) 
      collect_particles_point (point, p, multi_level_parts[l]);
    assign_particles(multi_level_parts[l], s);
  }
  multigrid_restriction(refl);
  return dep;
}

void cleanup_multi_level(Particles p) {
  free(multi_level_parts);
  for (scalar s in refl)
    free_scalar_data(s);
  delete (refl);
  free(refl);
  refl = NULL;

}
