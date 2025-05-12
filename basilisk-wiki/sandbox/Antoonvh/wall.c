/**
# A screen to block flow. 

Three experiments are performed,

1. A local surface heating in a stable quiescent atmoshphere  
2. 1 + a mean left-to-right flow  
3.  2 + a screen  

![Buoyancy and flow tracers for the three different experiments](wall/mov.mp4)
 */
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

scalar b[], * tracers = {b};
Particles parts;

b[bottom] = neumann(0.1*exp(-sq(x))/(mu.y[])); // A local heat source
b[top] = dirichlet (y);

face vector av[];
const face vector muc[] = {5e-3, 5e-3};

int maxlevel = 8, run_nr;

int main() {
  periodic (left);
  L0 = 20;
  X0 = -L0/2;
  mu = muc;
  a = av;
  N = 128;
  run();
  run_nr++;
  run();
  run_nr++;
  run();
}

double ux = 0.2;

event init (t = 0) {
  parts = init_tp_square (10, 0, 2, 3);
  foreach() {
    b[] = y;
    if (run_nr > 0)
      u.x[] = ux;
  }
}

event acceleration (i++) {
  foreach_face(y)
    av.y[] = (b[] + b[0,-1])/2.;
}

event tracer_diffusion (i++) 
  diffusion (b, dt, mu);
/**
## Damping layer and a fictious wall.
 */
#define WALL (y < 3 && fabs(x + 5) < 0.5)
event damp (i++) {
  foreach() {
    if (y > L0/2) 
      u.y[] *= exp(-(y - L0/2)/5.);
    if (run_nr == 2)
      if (WALL)
	foreach_dimension()
	  u.x[] *= exp(-dt*10);
  }
  boundary ((scalar*){u});
}

event mov (t += 0.5) {
  squares ("b", linear = true, min = -0.1, max = 4.1);
  scatter (parts);
  if (run_nr == 2) {
    scalar Wall[];
    foreach()
      Wall[] = WALL? 0 : nodata;
    squares ("Wall", map = gray, min = 0, max = 1);
  }
  mirror ({0,-1}) {
    cells();
  }
  save ("mov.mp4");
}

event adapt (i++) 
  adapt_wavelet ({b, u}, (double[]){0.1, 0.05, 0.05}, maxlevel, 5);

event stop (t = 100);