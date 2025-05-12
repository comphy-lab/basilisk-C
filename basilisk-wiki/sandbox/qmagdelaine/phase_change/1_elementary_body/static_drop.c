/**
# Drop evaporation

We investigate the evaporation of a single droplet in a still and relatively
dry environment.

![Evaporating droplet with the vapour concentration
field](static_drop/video_static_drop.mp4)

We define the geometrical, temporal and resolution parameters: */

#define R0 1.
#define L 10. // size of the box

#define MIN_LEVEL 5
#define LEVEL 10
#define MAX_LEVEL 12
#define dR_refine (2.*L0/(1 << LEVEL))

#define F_ERR 1e-10

#define T_END 600.
#define DT_MAX 10.
#define DELTA_T 10. // for measurements and videos

/**
We use an advection solver and the functions defined in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h): */

#include "axi.h"
#include "../advection_Q.h"
#include "../elementary_body.h"
#include "view.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

/**
## Physical parameters

We non-dimensionalise the problem with the initial radius of the drop
$R$, the diffusion coefficient of liquid vapour in air $D$ and the
saturation concentration of water in air $c_s$ (this quantity having
the dimension of a density $M\cdot L^{-3}$).

Thus the diffusion equation for the vapour in the non-dimensional
space reads:
$$
\frac{\partial c}{\partial t} = \Delta c,
$$
appropriately completed with Dirichlet conditions at the surface of
the drop and at infinity:
$$
\left\{
\begin{array}{ll}
c = 1 & \text{at the drop surface,}\\
c \to c_\infty/c_s & \text{far from the drop.}
\end{array}
\right.
$$
As the interface recedes at a (dimensional) velocity $v_\text{e} =
\frac{D}{\rho}\nabla c \sim \frac{D}{\rho}\frac{c_0}{R} \equiv V$, a
Peclet number Pe$= \frac{VR}{D} = \frac{c_s}{\rho}$ also enters in the
problem description. Typically, for the problem of a water droplet
evaporating in dry air, this Peclet number is $O(10^{-5})$, meaning
that the problem is really dominated by diffusion. Here, we choose
density ratio equal to $10^{-3}$, and set the relative humidity of the
surroundings to 20\%.

We need time factor to set the Dirichlet condition, its role is specified in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h). */

#define vapor_peclet 1e-3
#define D_V 1.
#define vcs 1.
#define cinf 0.2
#define dirichlet_time_factor 10.

/**
We allocate several scalar fields to describe both the
interface and the concentration field. */

scalar f[], vapor[];
scalar * interfaces = {f}, * tracers = {vapor};

/**
Thanks to symmetry, we only solve a quarter of the domain, and
requires the concentration to drop at its asymptotic value *at
infinity* (that is, at the box boundary) */

vapor[right] = dirichlet(cinf);
vapor[top]   = dirichlet(cinf);

/**
The main function of the program, where we set the domain geometry to
be ten times larger than the drop: */

int main() {
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);
  DT = DT_MAX;
  run();
}

/**
The initial position of the interface is defined with this function: */

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

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
  foreach()
    vapor[] = f[]*vcs + (1. - f[])*cinf;
  foreach_face()
    uf.x[] = 0.;
  boundary({vapor, uf});
  CFL = 0.2;
}

/**
## Evaporation velocity

The velocity due to evaporation is computed in the *stability()* event to take
into account this velocity in the CFL condition. */

event stability (i++) {
  vapor.D = D_V;
  vapor.peclet = vapor_peclet;
  vapor.inverse = true;
  phase_change_velocity (f, vapor, uf);
  boundary((scalar*){uf});
}

/**
After the *vof()* event, the evaporation velocity has to be erased. */

event tracer_advection (i++) {
  foreach_face()
    uf.x[] = 0.;
  boundary((scalar*){uf});
}

/**
## Diffusion with immersed dirichlet condition

The concentration field diffuses at each timestep. We need for that the maximal
level in the simulation. */

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
  adapt_wavelet ({f, vapor}, (double[]){1e-3, 1e-3},
		             minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif

/**
## Post-processings and videos

We now juste have to write several post-processing events to save the shape of
the drop, its effective radius and the vapor concentration profile. */

event interface (t = 3.*T_END/4.) {
  static FILE * fpshape = fopen("shape", "w");
  output_facets (f, fpshape);
  fflush(fpshape);
}

event outputs (t = 0.; t += max(DELTA_T, DT); t <= T_END) { 

  double effective_radius;
  effective_radius = pow(3.*statsf(f).sum, 1./3.);

  fprintf (stderr, "%.17g %.17g\n", t, effective_radius);
  fflush(stderr);
  if (t <= T_END - 20.) {
    static FILE * fpvapor = fopen("vapor_profile", "w");
    static FILE * fpvaporresc = fopen("vapor_profile_resc", "w");
    for (double y = 0.; y <= 5.; y += 0.01)
      fprintf (fpvapor, "%g %g\n", y, interpolate (vapor, 0., y));
    for (double y = effective_radius; y <= 5.; y += 0.01)
      fprintf (fpvaporresc, "%g %g %g\n", effective_radius, y/effective_radius,
               interpolate (vapor, 0., y));
    fprintf (fpvapor, "\n");
    fprintf (fpvaporresc, "\n");
    fflush(fpvapor);
    fflush(fpvaporresc);
  }
  
  /**
  We create a video with the concentration in the vapor phase. */

  scalar vapor_draw[];
  foreach() {
    f[] = clamp(f[], 0., 1.);
    vapor_draw[] = - vapor[];
  }
  boundary({f, vapor_draw});

  view (fov = 15, width = 640, height = 640, samples = 1, relative = false,
        tx = 0., ty = 0., bg = {BG, BG, BG});
  clear();
  draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
           fc = {BG, BG, BG});
  squares ("vapor_draw", min = - vcs, max = vcs, linear = false,
           map = cool_warm);
  mirror (n = {1., 0., 0.}, alpha = 0.) {
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
             fc = {BG, BG, BG});
    squares ("vapor_draw", min = - vcs, max = vcs, linear = false,
             map = cool_warm);
  }
  mirror (n = {0., 1., 0.}, alpha = 0.) { // vapor
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
             fc = {BG, BG, BG});
    squares ("vapor_draw", min = - vcs, max = vcs, linear = false,
             map = cool_warm);
    mirror (n = {1., 0., 0.}, alpha = 0.) {
        draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
                 fc = {BG, BG, BG});
        squares ("vapor_draw", min = - vcs, max = vcs, linear = false,
                 map = cool_warm);
    }
  }
  save ("video_static_drop.mp4");
}

/**
# Results

## Interface shape

We first check is the shape of drop remains spherical.

~~~gnuplot Shape of the interface
set style line 1 pt 7 ps 0.7
darkgray="#666666"
blue="#5082DC"
turquoise="#008C7D"
orange="#FF780F"
raspberry="#FA0F50"

set size square
set xlabel "x"
set ylabel "y"
set object 1 circle at 0,0 size first 0.4726 fc rgb darkgray
plot 'shape' w l lc rgb blue t 'simulation'
~~~

~~~gnuplot Local radius in function of the angular position
unset object 1
set xlabel "angle"
set ylabel "radius"
plot 'shape' u atan(1, 2):(sqrt($1*$1 + $2*$2)) ls 1 lc rgb blue t 'local radius vs angular position'
~~~

## Radius evolution: the $\text{d}^2$ law

This problem was first investigated independently by Maxwell and
Langmuir. They found that quite counter-intuitively, the evaporation
mass rate does not scale with the surface of the drop, but with its
radius (as a consequence of evaporation being a diffusion-limited
process).

It follows that the squared-radius of the drop decreases linearly with
time (so-called "$\text{d}^2$ law"), according to:
$$
r^2(t) = r_0^2 - 2 \frac{D(c_0-c_\infty)}{\rho}t
$$

Here is represented the time evolution of the effective radius of the
drop, compared with the theoretical prediction:

~~~gnuplot Time evolution of the squared radius $r^2(t)$
unset size
set xlabel "t"
set ylabel "R²"
plot [0:590][0:1] 'log' u 1:($2*$2) ls 1 lc rgb blue t 'simulation', \
     1.-0.002*0.8*x lw 1.5 lc rgb orange t 'inifinite model'
~~~

## Concentration profiles

Now checking the concentration profiles vs time we have

~~~gnuplot Raw concentration profiles vs time
set xlabel "r"
set ylabel "c"
plot 'vapor_profile' w l lc rgb blue t 'simulation'
~~~

The theoretical concentration profiles follow:
$$
c(r,t) = \left(\frac{R(t)}{r}\right) (c_0-c_\infty) + c_\infty
$$
This suggests to replot the concentration profiles versus the rescaled
variable $r/R(t)$:

~~~gnuplot Tentative rescaled concentration profiles
set xlabel "r/R"
set ylabel "c"
plot 'vapor_profile_resc' u 2:3 w l lc rgb blue t 'simulation', \
     0.8*(1./x) + 0.2  lw 1.5 lc rgb orange t 'infinite model'
~~~

This is not entirely satisfying.

The reason is that our infinity is 10; confinement effects, if not
dominant, are nonetheless not entirely negligible. Working out these
effects, we get:
$$
c(r,t) = c_0 + \frac{c_0-c_\infty}{1-\frac{R(t)}{R_\infty}}
\left(\frac{R(t)}{r} - 1\right)
$$

suggesting a slightly different rescaling:

~~~gnuplot Better rescaling for the concentration profiles, accounting for confinement effects
set xlabel "r/R"
set ylabel "c_r"
plot 'vapor_profile_resc' u 2:(($3-1.)*(1-$1/10.)+0.8) w l lc rgb blue t 'simulation', \
     0.8*(1./x) lw 1.5 lc rgb raspberry t 'confined model'
~~~

But wait. If taking confinement effects into account yield more neat
agreement, we should modify our $\text{d}^2$ law accordingly? Let's
give it a try. The confined $\text{d}^2$ law reads:
$$
R^2(t) - \frac{2}{3R_\infty} R^3(t) = R_0^2 - \frac{2}{3R_\infty}
R_0^3 - 2 \frac{D(c_0-c_\infty)}{\rho}t
$$

~~~gnuplot A better $\text{d}^2$ law
set xlabel "t"
set ylabel "R² - 2/3 R³/L_0"
plot [0:590][0:1] 'log' u 1:($2*$2-2./3./10.*$2*$2*$2) w p pt 7 ps 0.5 lc rgb blue t 'simulation', \
     1.-2./30.-0.002*0.8*x lw 1.5 lc rgb raspberry t 'confined model'
~~~

That's it!

# References

~~~bib
@article{maxwell1878,
  title={Diffusion},
  author={Maxwell, James Clerk},
  journal={Encyclopedia britannica},
  volume={7},
  pages={214--221},
  year={1878}
}

@article{langmuir1918,
  title={The evaporation of small spheres},
  author={Langmuir, Irving},
  journal={Physical review},
  volume={12},
  number={5},
  pages={368},
  year={1918}
}
~~~
*/
