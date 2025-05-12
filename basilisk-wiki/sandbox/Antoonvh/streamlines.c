/**
# Draw stream lines

On this page we visualize a 2D flow field, computed with `centered.h`,
via the vorticty field, tracer particles and streamlines. The movies
visualize a flow in a lid-driven cavity.

![The vorticity field](streamlines/vorticity.mp4)

![The tracer particles](streamlines/particles.mp4)

![The stream lines](streamlines/lines.mp4)

Each variant looks good and helps visualize the flow.
*/

#include "navier-stokes/centered.h"
#include "view.h"
#define BVIEW 1
#include "particles.h" //in /sandbox/Antoonvh/

#define UTOP (-(x - X0)*(x - (X0 + L0)))

u.t[top]  = dirichlet(UTOP); // A lid driven cavity

int main() {
  const face vector muc[] = {1./100., 1./100.};
  mu = muc;
  X0 = Y0 = -L0/2;
  N = 32;
  run();
} 

event init (t = 0) {
  init_particles_in_cells();
  refine (level <= 5);
}
/**
## Movie maker event

We generate the three aforementioned visualizations; 
*/ 
scalar psi[]; // In global scope for faster convergence.
event movie (t += 0.1; t <= 20) {
  scalar omega[];
  vertex scalar stream[];
  // Vorticity:
  vorticity (u, omega);
  squares ("omega");
  save ("vorticity.mp4");
  // Particles
  scatter (loc);
  box();
  save ("particles.mp4");
  // Stream lines:
  psi[bottom] = dirichlet (0);
  psi[top]    = neumann (-UTOP);   //Also a consistent BC for psi.
  psi[left]   = dirichlet (0);
  psi[right]  = dirichlet (0);
  poisson (psi, omega);
  boundary ({psi});
  foreach_vertex()
    stream[] = (psi[0,-1] + psi[-1,-1]
	      + psi[] + psi[-1])/4; 
  isoline ("stream", n = 14);
  box();
  save ("lines.mp4"); 
}
