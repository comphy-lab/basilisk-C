/**
# Oscillations of an evaporating drop

We simulate the oscillations of an evaporating droplet in a dry environment.
Just as [falling_drop.c](falling_drop.c) and [blown_drops.c](blown_drops.c),
it comes as a complement of [static_drop.c](static_drop.c) to show how the
implementation of the evaporation written in
[elementary_body.h](../elementary_body.h) can be easily coupled with
[Navier-Stokes](/src/navier-stokes/centered.h). The chosen parameters are not
physical at all.

![Evaporating droplet with the x-velocity (upper-left), the y-velocity field
(lower-left), vapour concentration (upper-right) and the pressure (lower-right)](vibrating_drop/video_vapor_p_ux_uy.mp4)

We define the geometrical, temporal and resolution parameters. The drop will be
initialized as an ellipse. Its long and short radii will be $a_e\, R0$ and
$b_e\, R0$, respectively. */

#define L 10. // size of the box
#define R0 1. // average initial radius

#define a_e 1.5
#define b_e (sqrt(1/a_e))

#define MIN_LEVEL 5
#define LEVEL 7
#define MAX_LEVEL 8

#define F_ERR 1e-10

#define T_END 40.
#define DELTA_T (T_END/400.) // for measurements and videos

/**
We use a Navier-Stokes solver and the functions defined in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h): */

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
see [static_drop.c](/sandbox/qmagdelaine/phase_change/1_elementary_body/static_drop.c)
and [elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h). */

#define vapor_peclet 1e-2
#define D_V 1.
#define vcs 1.
#define cinf 0.
#define dirichlet_time_factor 10.

#define RHOR 100.
#define MUR 10.

/**
The VOF tracer describing the interface, $f$, is already defined by
[two-phase.h](/src/two-phase.h). We allocate a second scalar field to describe
vapor concentration field and place it in the list of tracers of
[tracer.h](/src/tracer.h). */

scalar vapor[];
scalar * tracers = {vapor};

/**
Thanks to symmetry, we only solve a quarter of the domain, and requires the
vapour concentration at its asymptotic value *at infinity* (that is, at the box
boundary): */

vapor[right] = dirichlet(cinf);
vapor[top]   = dirichlet(cinf);

/**
We set the condition on the velocity which will ensure the wind: */

/**
In the main function of the program, we set the domain geometry to
be ten times larger than the drop, the surface tension, the gravity, the
densities and the viscoties of both phases. */

int main() {
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);

	f.sigma = 1.;
  rho1 = 1.;
  rho2 = rho1/RHOR;
  mu1 = 1e-3;
  mu2 = mu1/MUR;

  run();
}

/**
The initial position of the interface is defined with an ellipse: */

#define circle(x, y) (sq(R0) - sq(x/a_e) - sq(y/b_e))

/**
Before the first step, we initialize the concentration field (after having
refined the grid around the future interface): $c_s$ in the drop and 
$c_\infty$ in the vapor. $\mathbf{u}_f$ is set to zero. */

event init (i = 0) {
  fraction (f, circle(x,y));
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
[elementary_body.h](../elementary_body.h) for
details. */

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
  adapt_wavelet ({f, vapor, u}, (double[]){1e-3, 1e-2, 1e-2, 1e-2},
		             minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif

/**
## Post-processings and videos

We now juste have to write several post-processing events to save tthe effective
radius of the drop. */

double previous_delta = R0*(b_e - 1.);
double previous_time = 0.;

event compute_radius (i++) {
  scalar y_position[];
  coord D = {0., 1.};
  position (f, y_position, D, add=false);
  
  double y_extend = 0.;
  foreach_boundary(left) {
    if (interfacial(point, f) && y_position[] != nodata)
      y_extend = max(y_extend, y_position[]);
  }
  
  double effective_radius = pow(3.*statsf(f).sum, 1./3.);
  
  double delta_radius = y_extend - effective_radius;
  double period = t - previous_time;
  double period_th = 2.35*pow(effective_radius, 1.5);
  if (sign(delta_radius) > sign(previous_delta) && period > 0.8*period_th) {
    previous_time = t;
    if (period < 1.5*period_th) {
      static FILE * fperiod = fopen("period", "w");
      fprintf (fperiod, "%g %g %g\n", t, effective_radius, period);
      fflush(fperiod);
    }
  }
  previous_delta = delta_radius;
  static FILE * fpstep = fopen("step_data", "w");
  fprintf (fpstep, "%.17g %.17g %.17g\n", t, effective_radius, y_extend);
}

event outputs (t = 0.; t += DELTA_T; t <= T_END) {

  double effective_radius = pow(3.*statsf(f).sum, 1./3.);

  scalar velocity_in_drop[];
  foreach()
    velocity_in_drop[] = f[]*u.x[];
  boundary({velocity_in_drop});
  double effective_velocity = statsf(velocity_in_drop).sum/statsf(f).sum;
  
  fprintf (stderr, "%.17g %.17g %.17g\n", t, effective_radius, effective_velocity);
  fflush(stderr);

  /**
  We create a video with the concentration in the vapor phase and the velocity
  in the x (vertical) direction. */

  scalar vapor_draw[];
  foreach() {
    f[] = clamp(f[], 0., 1.);
    vapor_draw[] = - vapor[];
  }
  boundary({f, vapor_draw});

  view (fov = 18., width = 640, height = 640, samples = 1, relative = false,
        tx = -0., ty = -0., bg = {BG, BG, BG});
  clear();
  draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
           fc = {BG, BG, BG});
  squares ("vapor_draw", min = - vcs, max = vcs, linear = false,
           map = cool_warm);
  mirror (n = {0, 1, 0}, alpha = 0.) {
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
    squares ("p", min = - 3., max = 3., linear = false, map = cool_warm);
  }
  mirror (n = {1, 0, 0}, alpha = 0.) {
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0,
             fc = {BG, BG, BG});
    squares ("u.x", min = - 1., max = 1., linear = false,
             map = cool_warm);
    mirror (n = {0, 1, 0}, alpha = 0.) {
      draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0,
               fc = {BG, BG, BG});
      squares ("u.y", min = - 1., max = 1., linear = false,
               map = cool_warm);
    }
  }
  save ("video_vapor_p_ux_uy.mp4");
}

/**
# Results

## Radius evolution

We compare the radius evolution with the theortical one of a drop in a still
environnement confined in a box of size 20. The confined $\text{d}^2$ law reads:
$$
R^2(t) - \frac{2}{3R_\infty} R^3(t) = R_0^2 - \frac{2}{3R_\infty}
R_0^3 - 2 \frac{D(c_0-c_\infty)}{\rho}t
$$

~~~gnuplot Radius
darkgray="#666666"
purple="#9E5DB6"
blue="#5082DC"
turquoise="#008C7D"
forest="#149632"
orange="#FF780F"
raspberry="#FA0F50"
set style line 1 pt 7 ps 0.7 lw 0.5

set terminal @PNG enhanced size 640,640 font ",8"
set output '_radius.png'

set xlabel "t"
set ylabel "R² - 2/3 R³/L_0"

plot 'log' u 1:($2**2-2./3./10.*$2**3) w p ls 1 lc rgb blue t 'simulation', \
     1.-2./30.-0.02*x lw 1.5 lc rgb raspberry t 'model for a confined and still drop'
~~~

The agreement is rough. I see three possibles explanation to this:
low resolution, high Péclet (not quasi-stationary anymore) or vibrations. 

~~~gnuplot Oscillations
set terminal @PNG enhanced size 640,640 font ",8"
set output '_oscillations.png'

set xlabel "t"
set ylabel "R_y(t) - R_{eff}(t)"

plot 'step_data' u 1:($3/$2-1.) w lp ls 1 ps 0.2 lc rgb blue t 'wobulation'
~~~

We see nice oscillations of the apex of the drop on the y-axis.

~~~gnuplot Oscillation period
set terminal @PNG enhanced size 640,640 font ",8"
set output '_periods.png'

set xlabel "R_{eff}(t)"
set ylabel "T"

plot 'period' u 2:3 w p ls 1 lc rgb turquoise t 'period', \
     2.35*sqrt(x*x*x) lw 1.5 lc rgb raspberry t 'scaling law'
~~~

At each time, the oscillation period is well predicted by the scaling law
$\displaystyle T \sim \sqrt{\frac{\rho\, R^3}{\sigma}}$.
*/
