/**
# A wave-wall collision experiment

It turns out, initializing a compact solitary wave is not easy.

Therfore, we take inspiration from the lab setup: Waves are generated
by a moving piston inside a watertank of depth $h$. With the help of
the acceleration of gravity ($g$) waves will form and propagate to
collide with the basin's walls.
 */
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
						      /**
No-slip side walls and set and the moving piston parameters are chosen
adhoc.
							      */
   
u.t[right] = dirichlet (0.);
u.t[left] = dirichlet (0.);

double U_max = 1.5, //Maximum velocity
  tau = 0.25,       //Time scale for moving
  xp = 2.5,         //Start location
  w = 0.25,         //Width
  h = 3.,           //Max. height
  hb = 1;           //Height above bottom floor
#define U_X (1.5*exp(-sq((t - 2)/tau)))
#define PISTON (w - (fabs(x - xp)) - (y > h) - (y < hb))
scalar pstn[];
/**
Furthermore, domain size and fluid properties are also chosen adhoc.
 */
int main() {
  L0 = 10;
  mu1 = 0.01;
  mu2 = 0.04;
  f.sigma = 0.01;
  rho1 = 1000.;
  rho2 = 1.;
  run();
}
/**
## initialization

Ideally, a proper wave should be initialized here to avoid the
overhead of the piston wave making process. Now virtually all effort
is targeted towards the piston and not towards the actual wave - basin
wall collision.
 */
event init (t = 0) {
  pstn.prolongation = pstn.refine = fraction_refine;
  fraction (f, 2. - y); //Water depth $h = 2$
}
/**
The force of gravity is included with $g = 10$
 */
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] = -10;
}
/**
The moving piston is implemented via Stephane's trick. Note that this
piston is leaky.
 */
event piston (i++; t < 10) {
  xp += dt*U_X;
  fraction (pstn, PISTON);
  foreach() {
    u.y[] *= (1 - pstn[]);
    u.x[] = u.x[]*(1 - pstn[]) + pstn[]*U_X;
  }
}
/**
The grid is adapted to accuracy represent the piston fraction field,
the water fraction field and the velocity components. 
 */
event adapt (i++)
  adapt_wavelet ((scalar *){pstn, f, u}, (double[]){0.1, 0.02, 0.05, 0.05} , 9);
/**

## Movie

   For systems using Bview, `movie.mp4` may be generated.

![The water, piston and grid structure (mirrored)](wave_wall/movie.mp4)

We can see vortex shedding from the piston, wave overtopping and no
significant jetting at the walls. Meaning that the setup could be much
improved.
 */

#if 1
#include "view.h"
event bview_movie (t += 0.1) {
  view (tx = -0.5);
  draw_vof ("pstn", filled = 1, fc = {0.2,0.2,0.2});
  draw_vof ("f", filled = 1, fc = {0.1,0.1,0.9});
  mirror ({0,-1})
    cells();
  save ("movie.mp4");
}
/**
Else, `f.mp4` will have to reaveal the dynamics.
 */
#else
event movie (t += 0.1) {
  foreach() {
    if (pstn[] > 0.5)
      pstn[] = -1;
  }
  output_ppm (f, file = "f.mp4", mask = pstn,
	      n = 512, box = {{0,0},{10,5}});
}
#endif
/**
At t = 10, the simulation is stopped.
 */
event stop (t = 10.);


