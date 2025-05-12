/**
# Falling evaporating drop

We simulate the evaporation of a single falling droplet in a relatively dry
environment. It comes as a complement of
[static_drop.c](/sandbox/qmagdelaine/phase_change/1_elementary_body/static_drop.c)
to show how the implementation
of the evaporation written in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h)
can be easily coupled with [Navier-Stokes](/src/navier-stokes/centered.h).
The chosen parameters are not physical at all.

![Evaporating droplet with the vertical velocity field (left) and the vapour
concentration field (right)](falling_drop/video_ux_vapor.mp4)

We define the geometrical, temporal and resolution parameters: */

#define R0 1.
#define L 10. // size of the box

#define MIN_LEVEL 6
#define LEVEL 7
#define MAX_LEVEL 9
#define dR_refine (2.*L0/(1 << LEVEL))

#define F_ERR 1e-10

#define T_END 10.
#define DELTA_T (T_END/100.) // for measurements and videos

/**
We use a Navier-Stokes solver and the functions defined in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h): */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tracer.h"
#include "./../elementary_body.h"
#include "view.h"
#include "reduced.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

/**
We define the physical parameters associated to the evaporation. For details
see [static_drop.c](/sandbox/qmagdelaine/phase_change/1_elementary_body/static_drop.c)
and [elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h). */

#define vapor_peclet 1.
#define D_V 0.01
#define vcs 1.
#define cinf 0.
#define dirichlet_time_factor 10.

/**
The VOF tracer describing the interface, $f$, is already defined by
[two-phase.h](/src/two-phase.h). We allocate a second scalar field to describe
vapor concentration field and place it in the list of tracers of
[tracer.h](/src/tracer.h). */

scalar vapor[];
scalar * tracers = {vapor};

/**
Thanks to symmetry, we only solve half of the domain, and requires the vapour concentration at its asymptotic value *at
infinity* (that is, at the box boundary): */

vapor[right] = dirichlet(cinf);
vapor[top]   = dirichlet(cinf);
vapor[left]  = dirichlet(cinf);

/**
In the main function of the program, we set the domain geometry to
be ten times larger than the drop, the surface tension, the gravity, the
densities and the viscoties of both phases. */

int main() {
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);

	f.sigma = 1e-2;
  rho1 = 1.;
  rho2 = rho1/500.;
  mu1 = 1./10.;
  mu2 = mu1/20.;

  G.x = -0.1;
  Z.x = 0.;

  run();
}

/**
The initial position of the interface is defined with this function: */

#define circle(x, y, R) (sq(R) - sq(x-7.) - sq(y))

/**
Before the first step, we initialize the concentration field (after having
refined the grid around the future interface): $c_s$ in the drop and 
$c_\infty$ in the vapor. $\mathbf{u}_f$ is set to zero. */

event init (i = 0) {
  #if TREE
    refine (level < MAX_LEVEL && circle(x, y, (R0 - dR_refine)) < 0.
            && circle(x, y, (R0 + dR_refine)) > 0.);
  #endif
  fraction (f, circle(x, y, R0));
  foreach() {
    u.x[] = 0.;
    vapor[] = f[]*vcs + (1. - f[])*cinf;
  }
  foreach_face()
    uf.x[] = 0.;
  boundary({vapor, u, uf});
}

/**
## Evaporation velocity

The velocity due to evaporation is computed in the *stability()* event to take
into account this velocity in the CFL condition. Before to modify $\mathbf{u}_f$
we save it in another face vector field, in order to recover it just after the
advection of the interface. */

face vector uf_save[];

event stability (i++) {

  face vector ev[];

  vapor.D = D_V;
  vapor.peclet = vapor_peclet;
  vapor.inverse = true;
  phase_change_velocity (f, vapor, ev);
  
  foreach_face() {
    uf_save.x[] = uf.x[];
    uf.x[] += ev.x[];
  }
  boundary((scalar*){uf});
}

/**
After the *vof()* event, the evaporation velocity has to be set back to the
*real* velocity. The evaporation velocity is not a flow velocity but just a 
displacement of the interface. */

event tracer_advection (i++) {
  foreach_face()
    uf.x[] = uf_save.x[];
  boundary((scalar*){uf});
}

/**
## Diffusion with immersed dirichlet condition

The concentration field diffuses at each timestep. To ensure nicely the
saturation of the interface, we need the minimal diffusive time of the
simulation, which the diffusion time across the smallest cell. Therefore, we
need the maximal level in the simulation, see
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h) for
details. */

event tracer_diffusion(i++) {

#if TREE
  int max_level = MAX_LEVEL;
#else
  int max_level = LEVEL;
#endif
  
  vapor.D = D_V;
  vapor.tr_eq = vcs;
  vapor.inverse = true;
  dirichlet_diffusion (vapor, f, max_level, dt, dirichlet_time_factor);
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f, vapor, u}, (double[]){1e-3, 1e-1, 1e-1, 1e-1},
		             minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif

/**
## Post-processings and videos

We now juste have to write several post-processing events to save tthe effective
radius of the drop. */

event interface (t = 3.*T_END/4.) {
  static FILE * fpshape = fopen("shape", "w");
  output_facets (f, fpshape);
  fflush(fpshape);
}

event outputs (t = 0.; t += DELTA_T; t <= T_END) { 

  double effective_radius = pow(3.*statsf(f).sum, 1./3.);
  fprintf (stderr, "%.17g %.17g\n", t, effective_radius);
  fflush(stderr);

  /**
  We create a video with the concentration in the vapor phase and the velocity
  in the x (vertical) direction. */

  scalar vapor_draw[], u_draw[];
  foreach() {
    f[] = clamp(f[], 0., 1.);
    vapor_draw[] = - vapor[];
    u_draw[] = - u.x[];
  }
  boundary({f, vapor_draw, u_draw});

  view (fov = 18., width = 640, height = 640,
        tx = 0., ty = -0.5, bg = {BG, BG, BG},
        quat = {0., 0., - sqrt(0.5), sqrt(0.5)});
  clear();
  draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0,
           fc = {BG, BG, BG});
  squares ("u_draw", min = - 1., max = 1., linear = true,
           map = cool_warm);
  mirror (n = {0, 1, 0}, alpha = 0.) {
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
             fc = {BG, BG, BG});
    squares ("vapor_draw", min = - vcs, max = vcs, linear = true,
             map = cool_warm);
  }
  save ("video_ux_vapor.mp4");
}
