/**
# Internal waves using an adaptive grid

Similar to the succesful simulation [using the
multigrid](internalwacesMG.c) we now want to use adaptive mesh
refinement for the simulation of internal waves.

## Numerical set-up

The set-up is very similar to the [multigrid
simulation](internalwacesMG.c). This run is altered so that it uses an
adaptive quadtree grid.
*/
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "view.h"  //For cells()
#include "navier-stokes/perfs.h"

scalar b[], * tracers = {b};
face vector av[];
double sqN = 1, omega;

b[top]    = neumann (sqN);
b[bottom] = neumann (-sqN);

int main() {
  omega = sqrt(1./2.);
  a = av;
  L0 = 30;
  X0 = Y0 = -L0/2;
  TOLERANCE = 1e-4;
  DT = 0.2/omega;
  p.prolongation = p.refine = refine_linear; //3rd-order interpolation (vertical)
  N = 256;
  run();
}

event init (t = 0) {
  foreach()
    b[] = sqN*y;
  boundary ({b});
}

event acceleration (i++) {
  coord del = {0, 1};
  foreach_face()
    av.x[] = del.x*((b[] + b[-1])/2 +
		   0.1*(sin(omega*t)*((sq(x) + sq(y)) < 1)));
}

event output (t += 0.5; t <= 75) {
  scalar grb[];
  foreach() {
    grb[] = 0;
    foreach_dimension()
      grb[] += sq((b[1] - b[-1])/(2*Delta));
    grb[] = sqrt(grb[]);
  }
  boundary ({grb});
  squares ("grb", linear = true);
  cells();
  save ("mov.mp4");
}
/**
## Adaptivity

Inspired by the equations of motion, we decide to refine upon the
discrete representation for the byoyancy and velocity-component
fields.
*/
event adapt (i++)
  adapt_wavelet ({b, u}, (double[]){0.01, 0.01, 0.01}, 8, 4);

/**
## Results
   
![Waves and grid](internalwacesAMR/mov.mp4)

~~~gnuplot Convergence
set xlabel 'iteration'
set ylabel 'mgp.i / mgp.nrelax'
plot 'perfs' u 3, '' u 4
~~~
*/
