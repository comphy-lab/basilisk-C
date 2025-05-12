/**
# Evaporation of a plane film in an open box

We investigate the evaporation of plane film in a still and relatively dry
environment.

We define the geometrical, temporal and resolution parameters: */

#define H0 1.
#define L 10. // size of the box

#define MIN_LEVEL 4
#define LEVEL 6
#define MAX_LEVEL 8
#define dH_refine (2.*L0/(1 << LEVEL))

#define F_ERR 1e-10

#define T_END 100.
#define DT_MAX 0.05
#define DELTA_T 1. // for measurements

/**
$\mathrm{d}t_{\mathrm{max}}$ is the maximal timestep. We have to choose the
largest which allow us to get a nice transiant state: 0.005 seems perfect, 0.05
seems really close, 0.1 is still great, 0.2 is ok, 0.5 is not correct anymore.

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
$h_0$, the diffusion coefficient of liquid vapour in air $D$ and the
saturation concentration of water in air $c_s$ (this quantity having
the dimension of a density $M\cdot L^{-3}$).

Thus the diffusion equation for the vapour in the non-dimensional
space reads:
$$
\frac{\partial c}{\partial t} = \Delta c,
$$
appropriately completed with Dirichlet conditions at the surface of
the liquid and at the top of the box:
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
Peclet number Pe$= \frac{VR}{D} = \frac{c_0}{\rho}$ also enters in the
problem description. Typically, for the problem of a water droplet
evaporating in dry air, this Peclet number is $O(10^{-5})$, meaning
that the problem is really dominated by diffusion. Here, we choose
density ratio equal to $10^{-3}$, and set the relative humidity of the
surroundings to 20\%.

We need time factor to set the Dirichlet condition, its role is specified in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h). */

#define vapor_peclet 1e-3
#define D_V 1.
#define vcs 1.0
#define cinf 0.2
#define dirichlet_time_factor 10.

/**
We allocate several scalars and vector fields to describe both the
interface and the concentration field. */

scalar f[], vapor[];
scalar * interfaces = {f}, * tracers = {vapor};

/**
There is no flux at the bottom, the left and the right: this requires neumann
boundary conditions for the concentration. We just have to specifie the
dirichlet condition at the top, where we assume that the vapore remain constant
at $c_\infty$: */

vapor[top] = dirichlet(cinf) ;

/**
The main function of the program, where we set the domain geometry to
be ten times larger than the film: */

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

#define plane(x, y, H) (H - y)

/**
Before the first step, we initialize the concentration field (after having
refined the grid around the future interface): $c_s$ in the film and 
$c_\infty$ in the vapor. $\mathbf{u}_f$ is set to zero. */

event init (i = 0) {
  #if TREE
    refine (level < MAX_LEVEL && plane(x, y, (H0 - dH_refine)) < 0.
            && plane(x, y, (H0 + dH_refine)) > 0.);
  #endif
  fraction (f, plane(x, y, H0));
  foreach() {
    vapor[] = (1. - f[])*cinf + f[]*vcs;
  }
  foreach_face()
    uf.x[] = 0.;
  boundary({vapor, uf});
  CFL = 0.2;
}

/**
## Evaporation velocity

The velocity due to evaporation is computed in the *stability()* event to take
into account this velocity in the CFL condition. */

double max_velocity;

event stability (i++) {

  /**
  The first step is used to make diffuse the vapor, because at the beginning
  of the evaporation process, the concentration gradient diverges. We
  artificially make the simulation begin at a time $t_0 < 0$. In order to do
  this, we define a first timestep very short (we remain at $t \sim 0$), during
  which the interface does not move but the vapor diffuses over a finite time
  larger than the timestep. We will take into account this time translation in
  the model. */

  if (i == 0)
    DT = 0.001; 
  else {
    DT = DT_MAX;
    vapor.D = D_V;
    vapor.peclet = vapor_peclet;
    vapor.inverse = true;
    phase_change_velocity (f, vapor, uf);
    boundary((scalar*){uf});
  }
  max_velocity = - statsf(uf.y).min;
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

  if (i == 0) {
    double dt_save = dt;
    dt = 0.05;
    dirichlet_diffusion (vapor, f, max_level, dt, dirichlet_time_factor);
    dt = dt_save;
  }
  else
    dirichlet_diffusion (vapor, f, max_level, dt, dirichlet_time_factor);
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f, vapor}, (double[]){1e-3, 1e-3},
		             minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif

/**
## Post-processings

We now juste have to write several post-processing events to save the average 
thickness of film, the evaporation velocity and the vapor concentration profile.
*/

event outputs (t = 0.; t += max(DELTA_T, DT); t <= T_END) {

  double mean_thickness = statsf(f).sum/L0;
  fprintf (stderr, "%.17g %.17g %.17g\n", t, mean_thickness, max_velocity);
  fflush(stderr);

  static FILE * fpvapor = fopen("vapor_profile", "w");
  for (double y = 0.; y <= L0; y += L0/1000.)
    fprintf (fpvapor, "%g %g %g\n", t, y, interpolate (vapor, 0., y));
  fprintf (fpvapor, "\n");
  fflush (fpvapor);
}

/**
# Results

## Concentration profiles

In this problem the diffusion is coupled to the movement of the liquid surface. We assume that the flux of vapor is only diffusive :
$$
j = - D\, \frac{\partial c}{\partial y}
$$
\para The mass conservation in the vapor phase leads to an equation of diffusion:
$$
\frac{\partial c}{\partial t} = - \frac{\partial j}{\partial y} \quad \text{then} \quad \frac{\partial c}{\partial t} = D\, \frac{\partial^2 c}{\partial y^2}
$$
The recession of the film surface compensates the diffusive flux:
$$
\rho\, \frac{\mathrm{d}h}{\mathrm{d}t} = - j(h, t) \quad \text{then} \quad	\rho\, \frac{\mathrm{d}h}{\mathrm{d}t} = \left.D\, \frac{\partial c}{\partial y}\right\vert_h
$$
These equations show that the evaporation time $\tau_e$ scales as $\frac{\tau}{\mu}$ where $\mu = \frac{\Delta c}{\rho} \sim Pe$, $\Delta c = c_s - c_0$ and $\tau \sim \frac{\ell^2}{D}$, the diffusion time over a distance $\ell$.
The movement of the surface is then much slower than the diffusion. A quasi-static approach allow us to decouple the diffusion and the recession of the liquid surface: we look for the permanent solutions of the equation of diffusion and find:
$$
c(y) = c_s - \frac{\Delta c}{L - h} \, (y - h)
$$ 
Here are represented the time evolution of the concentration profile:

~~~gnuplot Concentration profiles though time
set style line 1 pt 7 ps 0.7
darkgray="#666666"
blue="#5082DC"
turquoise="#008C7D"
orange="#FF780F"
raspberry="#FA0F50"
set xlabel "y"
set ylabel "c"
plot 'vapor_profile' u 2:3 w l lc rgb blue t 'simulation'
~~~

We converge quicly toward a linear profile. 

## Thickness evolution

The thickness can be computed:
$$
\frac{\mathrm{d}h}{\mathrm{d}t} = \left.\frac{D}{\rho}\, \frac{\partial c}{\partial y}\right\vert_h \approx - \frac{\Delta c}{\rho}\, \frac{D}{L} \quad \text{thus} \quad h(t) \approx h_0 - \frac{\Delta c}{\rho}\, \frac{D}{L}\, t
$$
Here is represented the time evolution of the surface velocity, compared with the theoretical prediction:

~~~gnuplot Time evolution of the interface velocity
set xlabel "t"
set ylabel "v"
plot [0:60] 'log' u 1:3 ls 1 lc rgb blue t 'simulation', \
     0.001*0.8/9. lw 1.5 lc rgb orange t 'open box model'
~~~

And here is represented the time evolution of the mean thickness of the film, compared with the theoretical prediction:

~~~gnuplot Time evolution of the film thickness $h(t)$
set xlabel "t"
set ylabel "h"
plot 'log' u 1:2 w p pt 7 ps 0.6 lc rgb blue t 'simulation', \
     1.-0.001*0.8*x/9. lw 1.5 lc rgb orange t 'open box model'
~~~

After a quick transient regime, the velity reach the predicted constant velocity. Because of this transient regime, the theoretical prediction is shifted.

## Description of the transient regime

The transient regime lasts the time the linear concentration profile takes place. Before it does we assume that the boundary condition is at infinity. In this case, the complete equation can be solved without futher approximation with the similitude variable:
$$
x = \frac{y - h(t)}{\sqrt{D\, t}}
$$
We finally find:
$$
h_\infty(t) \underset{0}{\approx} h_0 - 2\, \mu\, \sqrt{\frac{D\, t}{\pi}} \quad \text{and} \quad v_\infty(t) \underset{0}{\approx} - \mu\, \sqrt{\frac{D}{\pi\, t}}
$$
This regime is supposed to last the time the linear profile takes place, i.e. a time of diffusion $\tau$ which should scales as $\frac{L^2}{D}$. It seems that the numerical factor is very closed to one fourth. The previous model can be corrected:
$$
h(t) \approx h_\infty(\tau) - \frac{\Delta c}{\rho}\, \frac{D}{L}\, (t - \tau)
$$
Here is represented the time evolution of the surface velocity, compared with the new theoretical predictions:

~~~gnuplot Time evolution of the evaporation velocity $v(t)$.
set samples 10000
set xlabel "t"
set ylabel "v"
plot [0:60] [0:0.00035]'log' u 1:3 w p pt 7 ps 0.6 lc rgb blue t 'simulation', \
     0.001*0.8/9. lw 1.5 lc rgb orange t 'open box model', \
     0.001*0.8* sqrt(1./(pi * (x+0.05))) lw 1.5 lc rgb raspberry t 'infinite model'
~~~

And here is represented the time evolution of the mean thickness of the film, compared with the new theoretical predictions:

~~~gnuplot Time evolution of the film thickness $h(t)$.
set xlabel "t"
set ylabel "h"
plot [0:100]'log' u 1:2 w p pt 7 ps 0.6 lc rgb blue t 'simulation', \
     1. - 2.*0.8*0.001*sqrt(81./(4.*pi)) - 0.001*0.8/9. * (x - 81./4.) lw 1.5 lc rgb orange t 'open box model', \
     1. +  2.*0.8*0.001*(sqrt(0.05/pi) - sqrt((x + 0.05)/pi)) lc rgb raspberry lw 1.5 t 'infinite model'
~~~

*/
