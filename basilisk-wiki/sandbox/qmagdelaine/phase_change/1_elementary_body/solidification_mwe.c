/**
# Solidification with a plane interface

We investigate the solidification of a block of water, between two plates, one
at -20°C and another at +20°C. It serves at minimum working example to show
how the functions defined in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h)
can be used to handle temperature fields and solidification.

![Solification of a block of water with the temperature
field](solidification_mwe/video_solidification.mp4)

We define the geometrical, temporal and resolution parameters: */

#define H0 1.
#define L 10. // size of the box

#define MIN_LEVEL 3
#define LEVEL 5
#define MAX_LEVEL 7
#define dH_refine (2.*L0/(1 << LEVEL))

#define F_ERR 1e-10

#define T_END 1000.
#define DT_MAX 10.
#define DELTA_T 10. // for videos and measurements

/**
We use an advection solver and the functions defined in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h): */

#include "../advection_Q.h"
#include "../elementary_body.h"
#include "view.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

/**
# Physical parameters

We non-dimensionalise the problem with the initial thickness of the film
$h_0$, the diffusion coefficient of temperature in the ice $D_S$ and the
command temperature $T_0$. A diffusive timescale appears: $\tau_{D_S} =
\frac{h_0^2}{D_S}$.

Thus the diffusion equation for the temperature in the non-dimensional
space reads:
$$
\frac{\partial \theta}{\partial t} = \Delta \theta,
$$
appropriately completed with Dirichlet conditions for the non-dimensionalise temperature $\theta$, at the top and bottom:
$$
\left\{
\begin{array}{ll}
\theta = 1 & \text{at the top}\\
\theta = -1 & \text{at the bottom}\\
\theta = 0 & \text{at the interface}
\end{array}
\right.
$$

The interface recedes at a velocity (expressed with dimensional quantities)
$$
v_\text{I} = \frac{1}{\rho_S\, L_H}\,\left(\lambda_S\, \nabla T\vert_S
+ \lambda_L\, \nabla T\vert_L\right)
\sim \frac{\lambda}{\rho_S\, L_H}\frac{T_0}{h_0} \equiv V
$$
a Peclet number Pe$= \frac{V\, h_0}{D_S} =
\frac{\lambda\, T_0}{L_H\, \rho_S\, D_S} = \frac{c_m^S\, T_0}{L_H}$
enters in the problem description, just as for the evaporation. This ratio is around 0.25 for
liquid water and 0.12 for ice.

We need time factor to set the Dirichlet condition, its role is specified in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h). */

#define T_eq 0.
#define D_L 0.12
#define D_S 1.
#define TL_inf 1.
#define TS_inf -1.
#define peclet_L 0.253
#define peclet_S 0.124
#define dirichlet_time_factor 10.

/**
We allocate several scalar fields to describe both the
interface and temperature fields. */

scalar f[], temperature_L[], temperature_S[];
scalar * interfaces = {f}, * tracers = {temperature_L, temperature_S};

/**
There is no flux at the left and the right: this requires neumann
boundary conditions for the temperatures. We just have to specifie the
dirichlet condition at the top and bottom: */

temperature_L[top] = dirichlet(TL_inf) ;
temperature_S[bottom] = dirichlet(TS_inf) ;

/**
The main function of the program, where we set the domain geometry to
be ten times larger than the initial thickness of ice: */

int main()
{
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);
  DT = DT_MAX;
  run();
}

/**
The initial position of the interface is defined with this function: */

#define plane(x, y, H) (y - H)

/**
Before the first step, we initialize the temperature fields (after having
refined the grid around the future interface). */

event init (i = 0) {
  #if TREE
    refine (level < MAX_LEVEL && plane(x, y, (H0 - dH_refine)) > 0.
            && plane(x, y, (H0 + dH_refine)) < 0.);
  #endif
  fraction (f, plane(x, y, H0));
  foreach() {
    temperature_L[] = f[]*TL_inf;
    temperature_S[] = (1. - f[])*TS_inf;
  }
  foreach_face()
    uf.x[] = 0.;
  boundary({temperature_L, temperature_S, uf});
  CFL = 0.2;
}

/**
## Solidification velocity

The velocity due to solidification is computed in the *stability()* event to
take into account this velocity in the CFL condition. */

event stability (i++) {

  /**
  The first step is used to make diffuse the temperature, because at the
  beginning, the concentration gradient diverges at the interface. We
  artificially make the simulation begin at a time $t_0 < 0$. In order to do
  this, we define a first timestep very short (we remain at $t \sim 0$), during
  which the interface does not move but the vapor diffuses over a finite time
  larger than the timestep. */

  if (i == 0)
    DT = 0.001; 
  else {
  
    /**
    In the evaporation examples
    ([static_drop.c](/sandbox/qmagdelaine/phase_change/1_elementary_body/static_drop.c) and
    [static_film.c](/sandbox/qmagdelaine/phase_change/1_elementary_body/static_film.c)),
    the evaporation velocity was due to the diffusive flux of the vapor in the
    gaz phase. Here, the solification is an equilibrium between the heat
    transfers of each side of the interface. */
  
    DT = DT_MAX;
    temperature_L.D = D_L;
    temperature_L.inverse = false;
    temperature_L.peclet = peclet_L;
    temperature_S.D = D_S;
    temperature_S.inverse = true;
    temperature_S.peclet = peclet_S;
    for (scalar t in tracers) {    
      face vector tv[];
      phase_change_velocity (f, t, tv);
      foreach_face()
        uf.x[] += (t.inverse ? tv.x[] : - tv.x[]);
    }
    boundary((scalar*){uf});
  }
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

The temperature diffuses at each timestep. We need for that the maximal
level in the simulation. */

event tracer_diffusion(i++) {
  
  #if TREE
    int max_level = MAX_LEVEL;
  #else
    int max_level = LEVEL;
  #endif

  temperature_L.D = D_L;
  temperature_L.inverse = false;
  temperature_S.D = D_S;
  temperature_S.inverse = true;
  for (scalar t in tracers)  
    t.tr_eq = T_eq;

  if (i == 0) {
    double dt_save = dt;
    dt = 0.05;
    for (scalar t in tracers)  
      dirichlet_diffusion (t, f, max_level, dt, dirichlet_time_factor);
    dt = dt_save;
  }
  else {
    for (scalar t in tracers)  
      dirichlet_diffusion (t, f, max_level, dt, dirichlet_time_factor);
  }
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f, temperature_L, temperature_S}, (double[]){1e-3, 1e-3, 1e-3},
		             minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif

/**
## Post-processings and video

We now juste have to write several post-processing events to save the average 
thickness of the ice and the temperature profile.
*/

event outputs (t = 0.; t += max(DELTA_T, DT); t <= T_END) {

  double mean_thickness = L0 - statsf(f).sum/L0;
  fprintf (stderr, "%.17g %.17g\n", t, mean_thickness);
  fflush(stderr);

  scalar temperature[];
  foreach()
    temperature[] = f[]*temperature_L[] + (1. - f[])*temperature_S[];
  boundary ({temperature});

  static FILE * fptemperature = fopen("temperature_profile", "w");
  for (double y = 0.; y <= L0; y += L0/1000.)
    fprintf (fptemperature, "%g %g %g\n", t, y, interpolate (temperature, 0., y));
  fprintf (fptemperature, "\n");
  fflush (fptemperature);
  
  /**
  We create a video with the temperature. */
  
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});
  
  view (fov = 18, width = 640, height = 640, samples = 1, relative = false,
        tx = -0.5, ty = -0.5, bg = {BG, BG, BG});
  clear();
  draw_vof ("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 0);
  squares ("temperature", min = TS_inf, max = TL_inf, linear = false,
           map = cool_warm);

  save ("video_solidification.mp4");
}

/**
# Results

## Concentration profiles

Here are represented the time evolution of the concentration profile:

~~~gnuplot Temperature profiles though time
set style line 1 pt 7 ps 0.7
darkgray="#666666"
blue="#5082DC"
turquoise="#008C7D"
orange="#FF780F"
raspberry="#FA0F50"
set xlabel "y"
set ylabel "T"
plot [0:10][-1:1]'temperature_profile' u 2:3 w l lc rgb blue t 'simulation'
~~~

## Thickness evolution

The velocity of the front is
$$
v_\text{e} = \frac{1}{L_H}\,\left(\frac{\lambda_S}{\rho_S}\, \nabla \theta\vert_S
+ \frac{\lambda_L}{\rho_L}\, \nabla \theta\vert_L\right)
$$
At equilibrium the diffusion solution in each face are linear profiles:
$$
\theta_S = - \frac{1}{h} \quad \text{and} \quad \theta_L = \frac{1}{L_0 - h}
$$
and the velocity is null, which gives the thickness at equilibrium:
$$
h_{eq} = \frac{L_0}{1 + \alpha} \quad \text{with} \quad
\alpha = \frac{\lambda_L\, \rho_S}{\lambda_S\, \rho_L}
$$

Here is represented the time evolution of the mean thickness of ice, compared
with the theoretical prediction:

~~~gnuplot Time evolution of the ice thickness $h(t)$
set xlabel "t"
set ylabel "h"
plot [:][0:10]'log' u 1:2 w p pt 7 ps 0.6 lc rgb blue t 'simulation', \
     10./(1. + 0.241) lw 1.5 lc rgb raspberry t 'equilibrium solution'
~~~

*/
