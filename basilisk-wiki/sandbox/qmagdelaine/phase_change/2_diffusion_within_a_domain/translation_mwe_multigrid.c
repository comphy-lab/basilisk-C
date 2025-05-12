/**
# Diffusion in a drop in translation

We make a tracer diffuse in a drop which is a translation. Thereafter, **this
diffusive tracer will be called solute**.
This minimum working example illustrates the combined use of an
immersed no flux condition and a simple advection, without adaptivity.
It complements the mwe
[translation_mwe](/sandbox/qmagdelaine/phase_change/2_diffusion_within_a_domain/translation_mwe.c)
where the adaptivity limits the conservation of the solute, which is not a
surprise since I did not write the proper prolongation and restriction functions.

Here the conservation of the solute is simply limited by the tolerance of the
Poisson solver. I did not expect the tolerance to have any influence since the
scheme should be conservative at the machine accuracy. Nevertheless, even
without any immersed boundary, whithout any advection, either with symetric or
periodic boundary condition, the solute is not perfectly conserved, see this mwe:
[diffusion_poisson](/sandbox/qmagdelaine/phase_change/2_diffusion_within_a_domain/diffusion_poisson.c).
Even if this loss is small, it limits the conservation of the solute in this
example.

![Diffusion in a moving droplet, concentration
field](translation_mwe_multigrid/video_translation.mp4)

Options: */

#define VIDEO 1

/**
We define the geometrical, temporal and resolution parameters: */

#define L 5.
#define R0 1.

#define MIN_LEVEL 5
#define LEVEL 9
#define MAX_LEVEL 9

#define errsolute 1e-7
#define F_ERR 1e-10

/**
The tolerance of the Poisson solver determines how well the solute is conserved,
so we set it very low: */

#define MY_TOLERANCE 1e-9

/**
We set the end time of the simulation and the time step. */

#define T_END 3. //3
#define DT_MAX 5e-2 // 5e-2 
#define DELTA_T 5e-2 // 5e-2 // for measurements and videos

/**
We use an advection solver and functions defined in
[mixtures.h](/sandbox/qmagdelaine/phase_change/mixtures.h) on a multigrid: */

#include "grid/multigrid.h"
#include "./../advection_Q.h"
#include "vof.h"
#include "./../mixtures.h"
#include "view.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

/**
We set diffusion coefficient of the solute and the velocity of the drop: */

#define D_L 0.1
#define VELOCITY 0.7

#ifndef AXI
  #define AXI 0
#endif

/**
We allocate several scalar fields to describe both the interface and the
concentration field. */

scalar f[], solute[];
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
We set the boundary condition of the uniform velocity field: */

uf.n[right] = dirichlet(fm.n[ghost]*VELOCITY);
uf.n[left] = dirichlet(fm.n[ghost]*VELOCITY);

/**
In the main function of the program, we set the domain geometry: */

int main() {
  size (L);
  origin (-L/2., 0.);
  TOLERANCE = MY_TOLERANCE;
  DT = DT_MAX;
  N = 1 << LEVEL;
  init_grid (N);
  run();
}

/**
The initial position of the interface is defined with this function: */

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

/**
Before the first step, we initialize the drop with a patch of solute at its
center, and the velocity field. */

double initial_solute_mass;

event init (i = 0) {
  fraction (f, circle(x + (L/2. - 1.5*R0), y, R0));
  fraction (solute, circle(x + (L/2. - 1.5*R0), y, 0.4*R0));
  foreach_face(x)
    uf.x[] = fm.x[]*VELOCITY;
  boundary({f, solute, uf});
  initial_solute_mass = (AXI ? 2.*pi : 2.)*statsf(solute).sum;
  CFL = 0.2;
}

/**
## Advection

The advection velocity is set in the *stability()* event to take into account
this velocity in the CFL condition. */

event stability (i++) {
  foreach_face(x)
    uf.x[] = fm.x[]*VELOCITY;
  boundary ((scalar *){uf});
}

/**
The solute is advected consistently with the VOF tracer by associating
the former to the latter, see [vof.h](/src/vof.h) for the details. */

static scalar * interfaces_save = NULL;

event vof (i++) {
  boundary ({f, solute, uf});
 
  /**
  The solute concentration field is associated to the VOF tracer $f$. */

  solute.inverse = false;
  solute.gradient = minmod2;
  f.tracers = {solute};

  /**
  We call here the *vof_advection* function. If not, the association of tracers
  to $f$ doesn't work. I don't know why, but I wrote a mwe here:
  [f_tracers_mwe.c](/sandbox/qmagdelaine/phase_change/3_non_solenoidal_advection/f_tracers_mwe.c).
  */

  vof_advection ({f}, i);

  boundary ({f, solute});

  /**
  We set the list of interfaces to NULL to prevent [vof.h](/src/vof.h) to do a
  second advection. */

  interfaces_save = interfaces;
  interfaces = NULL;
}

event tracer_advection (i++) {
  interfaces = interfaces_save;
}

/**
## Diffusion with immersed no flux condition

For the advection, we need the *quantity* field, $f\, c$, as explained in
[vof.h](/src/vof.h), but for the diffusion, we need the *concentration* field
$c$. Consequently, we divide the solute concentration field by $f$ before the
diffusion step and multiply it back by $f$ after the diffusion for the next
advection step. */

event tracer_diffusion (i++) {
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  foreach()
    solute[] = (f[] > F_ERR ? solute[]/f[] : 0.);
  boundary({solute});
  solute.D = D_L;
  solute.inverse = false;
  no_flux_diffusion (solute, f, dt);
  foreach()
    solute[] *= f[];
  boundary({solute});
}

/**
## Post-processings and videos

We now juste have to write several post-processing events to save the total
quantity of solute. */

event outputs(t = 0.; t += max(DELTA_T, DT); t <= T_END) { 

  double total_solute = (AXI ? 2.*pi : 2.)*statsf(solute).sum;

  fprintf (stderr, "%.17g %.17g\n", t, total_solute - initial_solute_mass);
  fflush(stderr);

  /**
  We create a video with the concentration in the liquid phase. */

  #if VIDEO
    view (fov = 18., width = 640, height = 640, samples = 1, relative = false,
          tx = 0., ty = 0., bg = {BG, BG, BG});
    clear();
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
    squares ("solute", min = -1., max = 1., linear = false, map = cool_warm);
    mirror (n = {0, 1, 0}, alpha = 0.) {
      draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
      squares ("solute", min = -1., max = 1., linear = false, map = cool_warm);
    }
    save ("video_translation.mp4");
  #endif 
}

/**
# Results

~~~gnuplot Gain or loss of solute
set style line 1 pt 7 ps 0.7
blue="#5082DC"
raspberry="#FA0F50"
set key right center
set format y '%.1e'
set xlabel "t"
set ylabel "error: gain or loss of matter"
plot 'log' u 1:2 ls 1 lc rgb raspberry t "gain or loss"
~~~

The loss is really low, but as in
[diffusion_poisson](/sandbox/qmagdelaine/phase_change/2_diffusion_within_a_domain/diffusion_poisson.c),
we have to decrease the tolerance to 1e-11 to obtain machine-like accuracy
($\sim 1e-14$). */