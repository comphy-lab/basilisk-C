/**
# Diffusion in a drop in translation

We make a tracer diffuse in a drop which is a translation. This minimum working
example illustrates the combined use of an immersed no flux condition and a
simple advection.

![Diffusion in a moving droplet, concentration
field](translation_mwe/video_translation.mp4)

Options: */

#define VIDEO 1
#define ADAPT 1
#define FACE_FRACTION_REFINE 0

/**
We define the geometrical, temporal and resolution parameters: */

#define L 5.
#define R0 1.

#define MIN_LEVEL 5
#define LEVEL 8
#define MAX_LEVEL 8

#define F_ERR 1e-6

#define T_END 3. //3
#define DT_MAX 5e-2 // 5e-2 
#define DELTA_T 5e-2 // 5e-2  for measurements and videos

/**
We use an advection solver and functions defined in
[mixtures.h](/sandbox/qmagdelaine/phase_change/mixtures.h): */

#include "axi.h"
#include "../advection_Q.h"
#include "vof.h"
#include "../mixtures.h"
#include "view.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

#define D_L 0.1
#define VELOCITY 0.7

#ifndef AXI
  #define AXI 0
#endif

/**
We allocate several scalar fields to describe both the interface and the
concentration field. */

scalar f[], ink[];
scalar * interfaces = {f}, * tracers = NULL;

/**
We add the diffusion coefficient as an attribute to the scalar fields. It is
defined in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h)
but this file is not included in this minimal working example. */

attribute {
  double D;
}

/**
to get a nice uniform velocity, we precise its limit conditions: */

uf.n[right] = dirichlet(fm.n[ghost]*VELOCITY);
uf.n[left] = dirichlet(fm.n[ghost]*VELOCITY);

/**
In the main function of the program, we set the domain geometry: */

int main() {
  size (L);
  origin (-L/2., 0.);
  DT = DT_MAX;
  N = 1 << LEVEL;
  init_grid (N);
  run();
}

/**
The initial position of the interface is defined with this function: */

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

/**
Before the first step, we initialize the drop with a patch of ink at its center,
and the velocity field.
*/

double initial_ink_mass;

event init (i = 0) {
  fraction (f, circle(x + (L/2. - 1.5*R0), y, R0));
  fraction (ink, circle(x + (L/2. - 1.5*R0), y, 0.6*R0));
  boundary({f, ink});
  initial_ink_mass = (AXI ? 2.*pi : 2.)*statsf(ink).sum;
  foreach_face(x)
    uf.x[] = fm.x[]*VELOCITY;
  boundary((scalar *){uf});
  CFL = 0.2;
}

/**
## Advection

The advection velocity is set in the *stability()* event to take into account
this velocity in the CFL condition. */


event stability (i++) {
  foreach_face(x)
    uf.x[] = fm.x[]*VELOCITY;
  boundary((scalar *){uf});
}

static scalar * interfaces_save = NULL;

event vof (i++) {
  boundary ({f, ink, uf});
 
  ink.inverse = false;
  ink.gradient = minmod2;
  f.tracers = {ink};

  vof_advection({f}, i);

  boundary({f, ink});

  interfaces_save = interfaces;
  interfaces = NULL;
}

event tracer_advection (i++) {
  interfaces = interfaces_save;
}

/**
## Diffusion with immersed no flux condition

For the advection, we need the *quantity* field, $f\, c$, as explained in
[vof.h](/src/vof.h), but for the diffusion, we need the *concentration* field.
*/

event tracer_diffusion(i++) {

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  foreach()
    ink[] = (f[] > F_ERR ? ink[]/f[] : 0.);
  boundary({ink});
  ink.D = D_L;
  ink.inverse = false;
  no_flux_diffusion (ink, f, dt);
  foreach()
    ink[] *= f[];
  boundary({ink});
}

#if TREE && ADAPT
  event adapt (i++) {
    boundary({ink});
    adapt_wavelet ({f, ink}, (double[]){1e-3, 1e-3},
		               minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
  }
#endif

/**
## Post-processings and videos

We now juste have to write several post-processing events to save the total
quantity of ink. */

event outputs(t = 0.; t += max(DELTA_T, DT); t <= T_END) { 

  double total_ink = (AXI ? 2.*pi : 2.)*statsf(ink).sum;

  fprintf (stderr, "%.17g %.17g\n", t, total_ink - initial_ink_mass);
  fflush(stderr);

  /**
  We create a video with the concentration in the vapor phase. */

  #if VIDEO
    view (fov = 18., width = 640, height = 640, bg = {BG, BG, BG});
    clear();
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
    squares ("ink", min = -1., max = 1., linear = false, map = cool_warm);
    mirror (n = {0, 1, 0}, alpha = 0.) {
      draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
      squares ("ink", min = -1., max = 1., linear = false, map = cool_warm);
    }
    save ("video_translation.mp4");
  #endif 
}

/**
# Results

~~~gnuplot Gain or loss of ink
set style line 1 pt 7 ps 0.7
blue="#5082DC"
raspberry="#FA0F50"
set key right center
set format y '%.1e'
set xlabel "t"
set ylabel "error: gain or loss of matter"
plot 'log' u 1:2 ls 1 lc rgb raspberry t "gain or loss"
~~~

*/
