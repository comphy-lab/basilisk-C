/**
# Cosmological dynamics

Consider an otherwise empty Newtonian universe with small but heavy
objects. The force of gravity seems important. For this purpose we use
[inertial-particles.h]() and the Poisson solver for the gravity
potential $\Phi$.
 */

scalar Phi[];
vector dPhi[];
 
#include "poisson.h"

mgstats mgP;
/**
   The inertial particles definitions must be tuned for this special
   purpose. First, particles each have a mass $m$. Note that we may
   loose compatibility with [tracer-particles.h]() is special care if
   not taken.
*/
#ifndef ADD_PART_MEM
#define ADD_PART_MEM coord u; double m; long unsigned int tag;
#endif
#include "inertial-particles.h"

Particles cosmos;

scalar * _automatics_ = NULL; //Not needed

event defaults (t = 0) {
  foreach_dimension() {
    Phi[left]  = dirichlet (0);
    Phi[right] = dirichlet (0);
  }
}
/**
The particles accelerate according to the gravity potential gradient.

$$\mathbf{a}_p = -\nabla \Phi.$$
 */

mgstats get_Phi (void) {
  particle_boundary (cosmos);
  scalar density[];
  foreach()
    density[] = 0;
  foreach_particle_in(cosmos)
    foreach_point (p().x,p().y,p().z)
      density[] += p().m;
  foreach()
    density[] /= dv();
  return poisson (Phi, density);
}

coord p_acc (PA_inp inp) {
  particle p = inp.p;
  coord a;
  double xp = p.x; double yp = p.y; double zp = p.z;  
  foreach_dimension()
    a.x = -interpolate (dPhi.x, xp, yp, zp); 
  return a;
}

event gravity (i++) {
  mgP = get_Phi();
  gradients ({Phi}, {dPhi});
}

/**
## Usage 

* [The big bang](bigbang.c)

*/
