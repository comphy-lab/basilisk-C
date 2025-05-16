/**
# Locate particles in a 3 point neighborhood

~~~gnuplot Coarse cell can "see" far away particles
set size ratio -1
set key outside
plot 'out' w l t 'cells',					\
'log' w l t  'centre-to-particle' ,\
'particles' w labels
~~~
 */
#include "particle_reference.h"
scalar p[];
Particles parts;
int np = 10;

int main() {
  srand (1);
  // setup grid
  L0 = 2.;
  X0 = Y0 = -L0/2.;
  init_grid (32);
  unrefine (x > 0);
  output_cells();

  // Allocate and initialize particles
  parts = new_particles (np);
  foreach_particle() {
    foreach_dimension()
      p().x = noise(); 
  }
  // Two particles share a cell for testing:
  pl[parts][np-1].x = pl[parts][0].x + 0.04;
  pl[parts][np-1].y = pl[parts][0].y + 0.01;
  
  assign_particles (parts, p);
  assign_particles (parts, p);   // Again for testing...
  // Print neighborhood data:
  foreach() {
    coord cc = {x, y};
    foreach_neighbor(1) {
      foreach_particle_point(p, point) {
	fprintf (stderr, "%g %g\n%g %g\n\n", cc.x, cc.y, x, y);
      }
    }
  }

  FILE * fp = fopen ("particles", "w");
  foreach_particle()
    fprintf (fp, "%g %g %d\n", x, y, _j_particle);
  fclose (fp);

  //clean_up
  free_p();
  free_scalar_data (p);
}

