/**
![Turbudity currents are driven by particles. Image via [Octavio E. Sequeiros et al. (2010)](https://ascelibrary.org/doi/10.1061/%28ASCE%29HY.1943-7900.0000200)](https://ascelibrary.org/cms/attachment/a5c1fe40-ea5f-4d00-a9d1-3a8694cd7137/7.gif)

# Vulcanic ash settling in water

We use 64000 spherical sand particles to model the descend of ganules
in water.

![You can "see" the two-way coupling in action](ash/locs.mp4)

The movie displays a surprisingly pleasing bview bug as is renders only part of the grid slice...
 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#define TWO_WAY 1
#include "stokes-particles.h"
#include "view.h"
#include "scatter2.h"

u.t[bottom] = dirichlet (0.);
u.r[bottom] = dirichlet (0.);

Particles sand;

#define MUZ 1e-3

int main() {
  const face vector muc[] = {MUZ, MUZ, MUZ}; //mu Water
  periodic (left);
  mu = muc;
  G.y = -9.81;
  L0 = 0.1;
  X0 = Y0 = Z0 = -L0/2; 
  N = 64;
  run();
}
/**
Initially, the flow is steady and seeded with particles near the top
of the domain.
 */

event init (t = 0) {
  const scalar rhof[] = 1000;
  rho = rhof;
  int pni = 40;
  if (pid() == 0) {
    sand = new_inertial_particles (cube(pni));
    int n = 0;
    for (int i = 0; i < pni; i++) 
      for (int j = 0; j < pni; j++)  
	for (int k = 0; k < pni; k++) { 
	  pl[sand][n].x = (double)i/(100*pni);
	  pl[sand][n].y = (double)j/(100*pni) + L0/4;
	  pl[sand][n].z = (double)k/(100*pni);
	  pl[sand][n].u2.x = 2000;                      // approx. rho Silicon
	  pl[sand][n].u2.y = 0.05e-3 + 0.01e-3*noise(); // fine sand-grain size
	  n++;
	}
  } else
    sand = new_inertial_particles (0);
  particle_boundary (sand);
}
/**
## Bottom boundary

Particles bounce at the bottom boundary. They should be able to rest
 on there.
*/
double damp = 2;
event bottom_boundary (i++) {
  foreach_particle_in(sand) {
    if (y < Y0) {
      p().y = Y0 + Y0 - y;
      foreach_dimension()
	p().u.y /= damp;
      p().u.y *= -1;  
    }
  }
}

event adapt (i++) 
  adapt_wavelet ((scalar*){u}, (double[]){0.01, 0.01, 0.01}, 7, 4);

event mov (t += .02) {
  view (theta = t/10. - 0.5);
  scatter (sand, pc = {1, 1, 1});
  box();
  translate (z = -L0/2) {
    cells();
  }
  save ("locs.mp4");
}

event end (t = 10);
