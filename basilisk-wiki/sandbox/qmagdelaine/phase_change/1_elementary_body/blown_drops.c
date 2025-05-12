/**
# Drops evaporating in the wind

We simulate the evaporation of a five droplets in the wind. Just as
[falling_drop.c](falling_drop.c) and [vibrating_drop.c](vibrating_drop.c),
it comes as a complement of [static_drop.c](static_drop.c) to show how the
implementation of the evaporation written in
[elementary_body.h](../elementary_body.h)
can be easily coupled with [Navier-Stokes](/src/navier-stokes/centered.h).
The chosen parameters are not physical at all.

![Evaporating droplets with the vorticity field (top) and the vapour
concentration field (bottom)](blown_drops/video_ux_vapor.mp4)

We define the geometrical, temporal and resolution parameters: */

#define L 10. // size of the box
#define R0 1.

#define RHOR 10.
#define MUR 10.

#define MIN_LEVEL 4
#define LEVEL 7
#define MAX_LEVEL 9

#define F_ERR 1e-10

#define T_END 25.
#define DELTA_T (T_END/200.) // for measurements and videos

/**
We use a Navier-Stokes solver and the functions defined in
[elementary_body.h](../elementary_body.h): */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tracer.h"
#include "../elementary_body.h"
#include "view.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

/**
We define the physical parameters associated to the evaporation. For details
see [static_drop.c](static_drop.c)
and [elementary_body.h](../elementary_body.h). */

#define vapor_peclet 5e-2
#define D_V 0.3
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
We set the condition on the velocity which will ensure the wind: */

u.n[left]  = dirichlet(1.);
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);

/**
In the main function of the program, we set the domain geometry to
be ten times larger than the drop, the surface tension, the gravity, the
densities and the viscoties of both phases. */

int main() {
  size (2.*L);
  origin (-L, 0.);
  N = 1 << LEVEL;
  init_grid (N);

	f.sigma = 0.1;
  rho1 = 1.;
  rho2 = rho1/RHOR;
  mu1 = 1e-3;
  mu2 = mu1/MUR;

  run();
}

/**
The initial position of the interface is defined with this function: */

#define circle(x, y) (sq(R0) - sq(x) - sq(y))

/**
Before the first step, we initialize the concentration field (after having
refined the grid around the future interface): $c_s$ in the drop and 
$c_\infty$ in the vapor. $\mathbf{u}_f$ is set to zero. */

event init (i = 0) {
  fraction (f, max(max(max(circle(x-8., y), circle(x-4., y)), max(circle(x, y), circle(x+4., y))), circle(x+8., y)));
  foreach() {
    u.x[] = 0.;
    vapor[] = f[]*vcs + (1. - f[])*cinf;
  }
  boundary({vapor, u});
  foreach_face()
    uf.x[] = 0.;
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
  foreach() {
    f[] = (f[] < F_ERR ? 0. : f[] > 1. - F_ERR ? 1. : f[]);
    f[] = clamp(f[], 0., 1.);
  }
  foreach_face()
    uf.x[] = uf_save.x[];
  boundary({f, uf});
}

/**
## Diffusion with immersed dirichlet condition

The concentration field diffuses at each timestep. To ensure nicely the
saturation of the interface, we need the minimal diffusive time of the
simulation, which the diffusion time across the smallest cell. Therefore, we
need the maximal level in the simulation, see
[elementary_body.h](../elementary_body.h) for details. */

event tracer_diffusion (i++) {
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
  adapt_wavelet ({f, vapor, u}, (double[]){1e-3, 1e-1, 1e-2, 1e-2},
                 minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif

/**
## Post-processings and videos

We now juste write a post-processing event to save a video. */

event outputs (t = 0.; t += DELTA_T; t <= T_END) {

  /**
  We create a video with the concentration in the vapor phase and the velocity
  in the x (vertical) direction. */
  scalar omega[];
  vorticity (u, omega);

  scalar vapor_draw[];
  foreach() {
    f[] = clamp(f[], 0., 1.);
    vapor_draw[] = - vapor[];
  }
  boundary({f, vapor_draw});

  view (fov = 19., width = 640, height = 640, samples = 1, relative = false,
        tx = 0., ty = 0., bg = {BG, BG, BG});
  clear();
  draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0,
           fc = {BG, BG, BG});
  squares ("omega", min = - 2., max = 2., linear = false,
           map = cool_warm);
  mirror (n = {0, 1, 0}, alpha = 0.) {
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
             fc = {BG, BG, BG});
    squares ("vapor_draw", min = - vcs, max = vcs, linear = false,
             map = cool_warm);
  }
  save ("video_ux_vapor.mp4");
}
