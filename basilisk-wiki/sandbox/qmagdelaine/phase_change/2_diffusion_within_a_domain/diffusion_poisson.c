/**
# Simple diffusion - losses

We investigate the conservation of a diffusive tracer. The diffusion equation is
solved thanks to the Poisson solver [poisson.h](/src/poisson.h). Since the later
uses fluxes between the cells, one can imagine that on
regular grids, with symetric or periodic boundary conditions, the tracer should
be perfectly conserved independently of the tolerance. It is not the case and
here is a minimum working example.

![Diffusive patch](diffusion_poisson/diffusive_patch.mp4)

Geometry and resolution: */

#define LEVEL 8

#define MY_TOLERANCE 1e-3
#define PERIODIC 1

#define T_END 2.
#define DELTA_T 5e-2

#define L 3.
#define R0 0.5

#include "run.h"
#include "fractions.h"
#include "poisson.h"
#include "view.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

/**
We allocate the diffusive tracer. */

scalar ink[];

/**
In the main function of the program, we set the domain geometry to
be ten times larger than the drop: */

int main() {
#if PERIODIC
  foreach_dimension()
    periodic(right);
#endif
  size (L);
  origin (-L/2., -L/2.);
  N = 1 << LEVEL;
  init_grid (N);
  TOLERANCE = MY_TOLERANCE;
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

double initial_ink_mass;

event init (i = 0) {
  fraction (ink, circle(x, y, R0));
  boundary({ink});
  initial_ink_mass = statsf(ink).sum;
}

event tracer_diffusion(i++) {

  /**
  We use directly [poisson.h](/src/poisson.h), but we follow the lines of
  [diffusion.h](/src/diffusion.h). */

  dt = dtnext (DELTA_T);
  scalar b[], lambda[];
  foreach() {
    b[] = -ink[]/dt;
    lambda[] = -1./dt;
  }
  face vector D[];
  foreach_face()
    D.x[] = 1.;
  boundary ({ink, b, lambda, D});
  poisson (ink, b, D, lambda);
}

/**
We now just have to save the losses and a nice video. */

event outputs (t = 0.; t += DELTA_T; t <= T_END) { 

  double total_ink = statsf(ink).sum;

  fprintf (stderr, "%.17g %.17g\n", t, total_ink/initial_ink_mass - 1.);
  fflush(stderr);

  view (fov = 18, width = 640, height = 640, samples = 1, relative = false,
        tx = 0., ty = 0., bg = {BG, BG, BG});
  clear();
  squares ("ink", min = -1., max = 1., linear = false, map = cool_warm);
  save ("diffusive_patch.mp4");
}

/**
# Results

~~~gnuplot ink losses
set format y '%.2e'
plot 'log' u 1:2 w p pt 7 ps 0.5 lc rgb "coral" t "created ink"
~~~

The error is around 1e-5 with the standard tolerance of 1e-3, with periodic or
symetric boundary conditions. With a toleance of 1e-11, the error vanishes
(around 1e-14). */