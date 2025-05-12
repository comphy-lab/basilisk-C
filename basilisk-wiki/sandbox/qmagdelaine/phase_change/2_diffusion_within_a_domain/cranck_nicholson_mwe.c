/**
# Diffusion in a drop

We make a tracer diffuse in a drop, applying a no flux condition at its surface.
This minimum working example illustrates the use of an immersed no flux
condition and we use it to compare a fully implicit scheme and the
Cranck-Nicholson scheme.

![Diffusion in a drop. Left: Cranck-Nicholson scheme, right: fully implicite
scheme.](cranck_nicholson_mwe/video_diffusion.mp4)

Options: */

#define VIDEO 1
#define ADAPT 1
#define FACE_FRACTION_REFINE 0

/**
We define the geometrical, temporal and resolution parameters: */

#define L 2.
#define R0 1.

#define MIN_LEVEL 5
#define LEVEL 8
#define MAX_LEVEL 8

#define F_ERR 1e-6

#define T_END 3. //3
#define DT_MAX 5e-2 // 5e-2 
#define DELTA_T 5e-2 // 5e-2  for measurements and videos

/**
We use the functions defined in 
[mixtures.h](/sandbox/qmagdelaine/phase_change/mixtures.h): */

#include "axi.h"
#include "run.h"
#include "fractions.h"
#include "../mixtures.h"
#include "view.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

#ifndef AXI
  #define AXI 0
#endif

/**
We define the diffusion coefficient: */

#define Diff 0.1

/**
We allocate several scalar fields to describe the interface and the
concentration fields: *ink* which will diffuse with a fully implicit scheme
and *ink_cranck* which will diffuse with Cranck Nicholson scheme. */

scalar f[], ink[], ink_cranck[];
scalar * liquid_tracers = {ink, ink_cranck};

/**
We add two attributes to the scalar fields: the diffusion coefficient $D$ and
a boolean to choose to which phase the tracer is associated. They are
defined in [vof.h](/src/vof.h) and
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h)
but these files are not included in this minimal working example. */

attribute {
  double D;
  double inverse;
}

/**
In the main function of the program, we set the domain geometry. */

int main() {
  size (L);
  origin (0., 0.);
  DT = DT_MAX;
  N = 1 << LEVEL;
  init_grid (N);
  run();
}

/**
The initial position of the interface is defined with this function: */

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

/**
We use the refine function of the fractions to prolongate and refine the VOF
tracer: */

event defaults (i = 0) {
  f.prolongation = f.refine = fraction_refine;
}

/**
Before the first step, we initialize the drop with a patch of ink at its center.
*/

double initial_ink_mass;

event init (i = 0) {
  fraction (f, circle(x, y, R0));
  for (scalar lt in liquid_tracers)
    fraction (lt, circle(x, y, 0.5*R0));
  boundary({f});
  boundary(liquid_tracers);
  initial_ink_mass = statsf(ink).sum;
  #if AXI
    initial_ink_mass *= 4.*pi;
  #endif
}

/**
## Diffusion with immersed no flux condition
*/

event tracer_diffusion(i++) {

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  dt = dtnext (DT_MAX);
  for (scalar lt in liquid_tracers) {
    lt.D = Diff;
    lt.inverse = false;
  }
  no_flux_diffusion (ink, f, dt);
  no_flux_diffusion_cranck (ink_cranck, f, dt);
}

/**
To ensure a good conservation of matter, we need to use the *quantity* field
$f\, c$ instead of the *concentration* field $c$, during the adaptation of the
grid. */

#if TREE && ADAPT
  event adapt (i++) {
    foreach() {
      for (scalar lt in liquid_tracers)
        lt[] *= f[];
    }
    boundary(liquid_tracers);
    adapt_wavelet ({f, ink, ink_cranck}, (double[]){1e-3, 1e-3, 1e-3},
		               minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
    foreach() {
      for (scalar lt in liquid_tracers)
        lt[] = (f[] > F_ERR ? lt[]/f[] : 0.);
    }
    boundary(liquid_tracers);
  }
#endif

/**
## Post-processings and videos

We now juste have to write several post-processing events to save the total
quantity of ink. */

event outputs (t = 0.; t += max(DELTA_T, DT); t <= T_END) { 

  attribute {
    double total;
  }
  
  for (scalar lt in liquid_tracers) {
    scalar quantity[];
    foreach()
      quantity[] = lt[]*f[];
    boundary({quantity});
    lt.total = (AXI ? 4.*pi : 1.)*statsf(quantity).sum;
  }

  fprintf (stderr, "%.17g %.17g %.17g\n", t, ink.total - initial_ink_mass,
           ink_cranck.total - initial_ink_mass);
  fflush(stderr);

  #if VIDEO
    view (fov = 30, width = 640, height = 640, bg = {BG, BG, BG});
    clear();
    // concentration ink
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
    squares ("ink", min = -1., max = 1., linear = false, map = cool_warm);
    mirror (n = {1., 0., 0.}, alpha = 0.) {
      draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
      squares ("ink", min = -1., max = 1., linear = false, map = cool_warm);
    }
    mirror (n = {0., 1., 0.}, alpha = 0.) {
      draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
      squares ("ink_cranck", min = -1., max = 1., linear = false, map = cool_warm);
      mirror (n = {1., 0., 0.}, alpha = 0.) {
        draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
        squares ("ink_cranck", min = -1., max = 1., linear = false, map = cool_warm);
      }
    }
    save ("video_diffusion.mp4");
  #endif 
}

/**
# Results

~~~gnuplot Gain or loss of ink
set style line 1 pt 7 ps 0.7
blue="#5082DC"
raspberry="#FA0F50"
set key right bottom
set format y '%.1e'
set xlabel "t"
set ylabel "error: gain or loss of matter"
plot 'log' u 1:2 ls 1 lc rgb blue t "fully implicit", \
     'log' u 1:3 ls 1 lc rgb raspberry t "Cranck-Nickolson"
~~~

~~~gnuplot Gain or loss of ink (fully implicit)
set xlabel "t"
set ylabel "error: gain or loss of matter"
plot 'log' u 1:2 ls 1 lc rgb blue t "fully implicit"
~~~

The Cranck-Nicholson scheme allows to reduce the loss of matter if a constant
grid is used, but with an adaptative one, the fully implicit one is much better.
Moreover, on both constant and adaptative grids, some glitches appear in the
concentration field with the Cranck-Nicholson scheme (see the video at the top
of the page).
*/
