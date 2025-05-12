/**
# Film Boiling

The film boiling configuration consists in a superheated
solid wall, covered by a thin vapor layer, which is
perturbed by the phase change phenomena releasing bubbles.
The simulation setup used here was inspired by [Welch and
Wilson, 2000](#welch2000volume). The domain is initialized
with a uniform constant temperature equal to the saturation
value. The bottom wall is maintained at a specific
superheating temperature, causing the vapor layer to
expand and leading to the departure of bubbles, which
rise the liquid column unil the breakup with the free
surface.

![Evolution of the temperature field and the gas-liquid interface](filmboiling/movie.mp4)(width="800" height="600")
*/

/**
## Phase change setup

We do not use the Stefan flow shifting procedure. The
volume fraction source term is divided by the density of
the gas phase. We compute the interface gradients just
in gas-phase, and we use the same velocity used for the
volume fraction also for the tracers in phase 2: *TG*.
Defining the variable *WELCH* you can use the setup
described in the [paper](#welch2000volume).
*/

#define BOILING_SETUP
#define SOLVE_GASONLY

/**
## Simulation setup

We use the centered solver with the evaporation source term.
The velocity potential method is adopted to obtain a
divergece-free velocity extension for the VOF advection.
The temperature-gradient phase change model is employed,
since just the temperature field has to be solved.
*/

#if POTENTIAL
# include "navier-stokes/centered-evaporation.h"
# include "navier-stokes/velocity-potential.h"
#else
# include "navier-stokes/velocity-jump.h"
#endif
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "evaporation.h"
#include "temperature-gradient.h"
#include "tag.h"
#include "view.h"

/**
The characteristic length of the problem is the most
dangerous Taylor wavelength, defined as:
$$
  \lambda_0 = 2\pi
  \left(
  \dfrac{3\sigma}{(\rho_l - \rho_g)g}
  \right)^{1/2}
$$

Then, we declare variables necessary for the temperature
gradient model, the maximum level of refinement, and the
temperature of the solid wall.
*/

double wavelength = 0.0786844;
int maxlevel = 8;

double lambda1, lambda2, cp1, cp2, dhev, TIntVal;
double TL0, TG0, Twall, Tsat, gas_vol0, width;

/**
We set the boundary conditions for velocity, pressure
and temperature at the bottom, which is considered as
as superheated solid wall. The velocity potential *ps*
boundary conditions are set in order to be coherent with
the pressure boundary conditions.
*/

#if POTENTIAL
u.n[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);
p[bottom] = neumann (0.);
ps[bottom] = neumann (0.);

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
ps[top] = dirichlet (0.);
#else
u1.n[bottom] = dirichlet (0.);
u1.t[bottom] = dirichlet (0.);
u2.n[bottom] = dirichlet (0.);
u2.t[bottom] = dirichlet (0.);
p[bottom] = neumann (0.);
ps[bottom] = neumann (0.);
pg[bottom] = neumann (0.);

u1.n[top] = neumann (0.);
u1.t[top] = neumann (0.);
u2.n[top] = neumann (0.);
u2.t[top] = neumann (0.);
p[top] = dirichlet (0.);
ps[top] = dirichlet (0.);
pg[top] = dirichlet (0.);
#endif

f[bottom] = dirichlet (0.);
TG[bottom] = dirichlet (Twall);
TL[bottom] = dirichlet (Twall);

int main (void) {
  /**
  We set the values of the material properties of the
  fluids. */

  rho1 = 200., rho2 = 5.;
#ifdef WELCH
  mu1 = 0.1, mu2 = 0.005;
#else
  mu1 = 0.01, mu2 = 0.005;
#endif
  lambda1 = 40., lambda2 = 1.;
  cp1 = 400., cp2 = 200.;
#ifdef WELCH
  dhev = 1.e4;
#else
  dhev = 1.e5;
#endif

  /**
  We set the initial and wall temperatures. */

  Tsat = 500., Twall = 506.;
  TL0 = Tsat, TG0 = Tsat, TIntVal = Tsat;

  /**
  We set the surface tension coefficient, and we
  apply the gravity force using the [reduced.h](/src/reduced.h)
  approach. */

  f.sigma = 0.1;
  G.y = -9.81;

  /**
  The dimensions of the problems are a function of the
  wavelength. In principle, the vapor bubbles, should
  be positioned in a squared pattern separated by a 
  distance equivalent to $\lambda_0$. */

  size (3.*wavelength);
  init_grid (1 << maxlevel);
  run();
}

/**
We initialize the volume fraction field and the temperature
field accordingly. The amount of gas in the domain at the
beginning of the simulation is stored.
*/

#define sin(x,y)(y-wavelength/128*(4+cos(2*pi*x/wavelength)))

event init (i = 0) {
  width = 0.5*wavelength;
  mask (x > width ? right : none);
  fraction (f, sin(x,y));
#ifndef WELCH
  foreach()
    f[] = (y > 2./3.*L0) ? 0. : f[];
#endif
  foreach() {
    TL[] = TL0*f[];
    TG[] = TG0*(1. - f[]);
    T[]  = TL[] + TG[];
  }

  foreach (reduction(+:gas_vol0))
    gas_vol0 += (1. - f[])*dv();
}

/**
We write in the stdout the simulation time, the amount of
gas in the domain, the Nusselt number from this simulation:
$$
  Nu = \dfrac{\lambda_0}{w \left(T_{wall}-T_{sat}\right)}
  \int_{0}^{w}
  \left.\dfrac{\partial T}{\partial y}\right\vert_{y=0} dw
$$
and the Nusselt number from the Berenson correlation:
$$
  Nu_B = 0.425
  \left(
  \dfrac{\rho_g(\rho_l-\rho_g)g\Delta h_{ev}}
  {\lambda_g \mu_g (T_{wall}-T_{sat})}
  \right)^{1/4} (\lambda_0)^{3/4}
$$
*/

event logfile (i++) {
  double gas_vol = 0.;
  foreach (reduction(+:gas_vol))
    gas_vol += (1. - f[])*dv();

  foreach()
    T[] = TL[] + TG[];

  double Nu = 0.;
  foreach_boundary (bottom, reduction(+:Nu)) {
    T[0,-1] = 2.*Twall - T[];
    Nu += (T[] - T[0,-1]);
  }
  Nu *= -wavelength/(Twall - Tsat)/width;

  double NuB = 0.425*pow((rho2*(rho1 - rho2)*9.81*dhev)
      /(lambda2*mu2*(Twall - Tsat)), 1./4.)*pow(wavelength, 3./4.);

  fprintf (stdout, "%f %f %f %f\n", t, gas_vol / gas_vol0, Nu, NuB);
  fflush (stdout);
}

/**
We remove small bubbles formed during the breakups that are
not important for the process under investigation, but that
give some problems during the solution of the diffusion
part of the temperature equation. */

event remove_droplets (i++) {
  remove_droplets (f, threshold=F_ERR, bubbles=true);
}

/**
We output a movie with the evolution of the temperature
field and the volume fraction facets. */

event movie (t += 0.01; t <= 6) {
  clear();
  view (ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = Tsat, max = Twall,
      map=blue_white_red, linear = true);
  mirror ({1,0}) {
    draw_vof ("f", lw = 1.5);
    squares ("T", min = Tsat, max = Twall,
        map=blue_white_red, linear = true);
  }
  save ("movie.mp4");
}

/**
## Results

We compare the evolution of the Nusselt number obtained
from the simulation with the Berenson correlation.
A displacement is expected since the correlation describes
a 3D system, while the simulation is 2D.

~~~gnuplot Evolution of the Nusselt number
set yr[0:60]

set xlabel "t [s]"
set ylabel "Nu [-]"

p "out" u 1:3 w l t "Results", "out" u 1:4 w l t "Berenson"
~~~

The amount of gas phase volume fraction increases according
to the phase change. The discontinuities correspond to the
bubbles release.

~~~gnuplot Evolution of gas volume fraction
reset
set xlabel "t [s]"
set ylabel "gas volume fraction [-]"

p "out" u 1:2 w l t "volume fraction"
~~~

## References

~~~bib
@article{welch2000volume,
  title={A volume of fluid based method for fluid flows with phase change},
  author={Welch, Samuel WJ and Wilson, John},
  journal={Journal of computational physics},
  volume={160},
  number={2},
  pages={662--682},
  year={2000},
  publisher={Elsevier}
}
~~~
*/
