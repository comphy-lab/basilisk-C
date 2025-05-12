/**
# Test for `pdump` and `prestore` 

Here we test if particles are dumped and restored correctly with
`_MPI`.

There will be a four particle lists in this test,
*/
#include "particle.h"
#include "view.h"
#include "scatter2.h"

Particles zeroth, first, second, third;

/**
   A helper function to view the particle locations and the `pid()` is
   defined below.
*/

void plot (char * fname, Particles * list) {
  foreach_P_in_list (list) {
    translate (z = 0.005)
      scatter (P, pc = {cos(P), sin(P), sin(P*1.4)});
    scatter (P, s = 50,
	     pc = {sin(pid()), cos(pid()), sin(pid()*2.4)});
  }
  translate (z = -0.01)
    cells();
  save (fname);
}

int main() {
  X0 = Y0 = -L0/2;
  init_grid (4);
  refine (x + y < 0 && level < 3);
  /**
     For testing purposes, the zeroth list is empty.
  */
  zeroth = new_particles (0);
  /**
     The first list contains particles at cell centers.
  */
  int n = 0;
  foreach()
    n++;
  first  = new_particles (n);
  place_in_cells (first);
  /**
     The second and third list are 25 and 30 particles in a square and
     circle, respectively.
  */
  if (pid() == 0) {
    second  = new_particles (25);
    third   = new_particles (30);
    place_in_square (second, (struct Init_P){5});
    place_in_circle (third,  (struct Init_P){30, 0.25, 0.25, L0/8.});
  } else {
    second  = new_particles (0);
    third   = new_particles (0);
  }
  Particles pall[4] = {zeroth, first, second, third};
  foreach_P_in_list (pall)
    particle_boundary (P);
  Particles Plist[3] = {zeroth, first, third};
  
  plot ("all.png",  pall);
  plot ("Plist.png", Plist);
  /**
     The particle locations with `list` center-colour coding and `pid()`
     ring colour.
     
     ![All particles](test_prestore/all.png)

     ![Particles in `Plist`](test_prestore/Plist.png)
     
     The particles in the `Plist` particles lists are dumped.
  */
  pdump ("dumplist", Plist);
  /**
     We create a new grid, randomize the old particles.
  */
  init_grid (4);
  refine (x + y > 0 && level < 3);
  foreach_particle() 
    foreach_dimension()
      p().x = 0.75*noise();
  plot ("pmess.png", pall);
  /**
     Particles are a mess before restoring

![Scrambled particles](test_prestore/pmess.png)

Now some data is restored:
   */
  
  prestore ("dumplist");
  /**
     The `Plist` may need be redefined for further usage.
   */
  // Plist[0] = zeroth = 0;
  // Plist[1] = first  = 1;
  Plist[2] = third  = 2;
  plot ("restored.png", Plist);
  /**
## The result: 

![Particles in `list`](test_prestore/restored.png)

Close inspection of the colour codes learns that everything went well!

   */
  free_p ();
}
