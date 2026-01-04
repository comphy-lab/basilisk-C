/**
![transition from laminar to turbulent](https://upload.wikimedia.org/wikipedia/commons/0/03/Laminar-turbulent_transition.jpg){width=300px}

# A rising plume

Herein we model the flow originating from a localized heat source. Such flows are often associated with smoke. As such, the flow is visualized with a volumetric rendering of the smoke concentration field.

![Volumetric rendering](heating/plume.mp4)

![grid refinement](heating/level.mp4)
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "utils.h"
#include "bwatch.h"

/**
The buoyancy and smoke fields are modelled seperately
*/

scalar b[], s[], * tracers = {b, s};
//size of heat source 
double siz = 0.3;

face vector av[];

double muv = 2e-4;
// Grid paramters
double be = 0.05, ue = 0.05;
int maxlevel = 10, minlevel = 5;


int main() {
  L0 = 100;
  X0 = Y0 = Z0 = -L0/2;
  a = av;
  DT = 0.1;
  N = 1 << minlevel;
  const face vector muc[] = {muv, muv, muv};
  mu = muc;
  run();
}

event init (t = 0) {
  for (scalar s in tracers)
    s.gradient = minmod2;
  refine ((sq(x) + sq(y + L0/3.) + sq(z)) < sq(20*siz) &&
	  level < maxlevel - 2);

  refine ((sq(x) + sq(y + L0/3.) + sq(z)) < sq(10*siz) &&
	  level < maxlevel - 1);
  refine ((sq(x) + sq(y + L0/3.) + sq(z)) < sq(2*siz) &&
	  level < maxlevel);
  foreach()
    for (scalar st in tracers)
      st[] = exp(-(sq(x) + sq(y + L0/3.) + sq(z))/sq(siz));

}

event acceleration(i++) {
  double tau = 1;
  foreach()
    for (scalar st in tracers)
      st[] += (1 - st[])*dt/tau*exp(-(sq(x) + sq(y + L0/3.) + sq(z))/sq(siz));
  
  foreach_face(y) 
    av.y[] = (b[] + b[0,-1])/2.;
}
/**
Only the buoyancy field is subjected to diffusion. 
*/
event tracer_diffusion (i++) {
  diffusion (b, dt, mu);
}

event adapt (i++) {
  adapt_wavelet  ({s, b, u}, (double[]){be/4., be, ue, ue, ue},
		  maxlevel, minlevel);
}

event mov (t += 0.5) {
  output_ppm (b, file = "b.mp4", n = 300);
  output_ppm (s, file = "s.mp4", n = 300);
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (lev, file = "level.mp4", n = 300,
	      min = minlevel - 1, max = maxlevel);
  
  static FILE * fp = popen ("ppm2mp4 plume.mp4", "w");
  watch (fov = 40, O = {1,1,100}, nx = 400, ny = 800);
  sphere (150, mat = {.dull = true, .col[0]  = 5, .col[1]  = 5, .col[2]  = 70,
		      .col2[0] = 5, .col2[1] = 5, .col2[2] = 5,
		      .n1 = (coord){0.1, 1, 0}});
  volume (s, sc = siz, mval = 1e-2, col = {255,255,255}, shading = 1); 
  store (fp);
}

event stop (t = 125) {
  dump ("smoke_dump");
}