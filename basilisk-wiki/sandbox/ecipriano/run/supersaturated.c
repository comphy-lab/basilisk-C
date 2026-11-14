/**
# Bubble Growth in a Supersaturated Solution

A spherical gas bubble is suspended in a liquid supersaturated with the
dissolved gas. The gas diffuses towards the interface and the bubble grows.
The setup is the one of [Gennari et al., 2022](#gennari2022), Sec. 4.4.1. The
bubble is pure gas and it is the *dispersed* phase, while the liquid is the
*continuous* phase, made of the dissolved gas and of water, the inert species
(see [evaporation.h](../src/evaporation.h)). The gas moves from the continuous
to the dispersed phase.

This is the mass transfer version of the [Scriven problem](scrivenproblem.c),
where a vapor bubble grows in a superheated liquid: the supersaturation of the
dissolved gas plays the role of the liquid superheating, the interface mass
fraction from Henry's law replaces the saturation temperature, and the problem
is solved by [evaporation.h](../src/evaporation.h) instead of
[boiling.h](../src/boiling.h).

The analytic solution was found by [Scriven, 1959](#scriven1959): the bubble
radius grows as
$$
  R(t) = 2\beta\sqrt{\mathcal{D}t}
$$
where the growth constant $\beta$ is given by
$$
  \Phi = 2\beta^3\exp(\beta^2 + 2\varepsilon\beta^2)
    \int_\beta^\infty x^{-2}\exp(-x^2 - 2\varepsilon\beta^3 x^{-1})\,dx,
  \quad
  \Phi = \dfrac{\rho_l}{\rho_g}\dfrac{\omega_{bulk} - \hat\omega}
    {1 - \hat\omega},
  \quad
  \varepsilon = \dfrac{\rho_l - \rho_g}{\rho_l}
$$
and $\hat\omega$ is the interface mass fraction of the dissolved gas, given by
Henry's law $\hat{c}_l = \hat{c}_g/He$, i.e. $\hat\omega = \rho_g/(He\rho_l)$
for a pure bubble. The bulk mass fraction is $\omega_{bulk} = \zeta\hat\omega$,
with $\zeta$ the saturation ratio. The concentration field is self-similar
$$
  \omega(r,t) = \omega_{bulk} - (\omega_{bulk} - \hat\omega)
    \dfrac{F(r/2\sqrt{\mathcal{D}t})}{F(\beta)},
  \quad
  F(s) = \int_s^\infty x^{-2}\exp(-x^2 - 2\varepsilon\beta^3 x^{-1})\,dx
$$
Since the analytic solution starts from a bubble of zero size, the simulation
starts at the time $t^*$ when the radius is $R_0$, with the concentration field
initialized from the analytic solution.

We run three levels of refinement up to $t = 4$ s, when the radius grows
from 5 to 9.1 mm, and we measure the convergence of the radius error. */

#include "axi.h"
#include "navier-stokes/low-mach.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"

#include <gsl/gsl_integration.h>
#pragma autolink -lgsl -lgslcblas

/**
### Boundary conditions

The bubble is centered in the bottom-left corner, exploiting the spherical
symmetry. Outflow boundary conditions are imposed on the top and right sides.
The boundary conditions of the species are set in *init*, since the species
fields are created in *defaults*. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

/**
### Problem data

The saturation ratio $\zeta$ can be changed at compilation time (2, 3, 4 in
the paper). */

#ifndef ZETA
# define ZETA 2.
#endif

#ifndef TMAX
# define TMAX 4.
#endif

int maxlevel, minlevel = 4;
double R0 = 0.005, He = 5., Dl = 2e-5;
double wint, wbulk, epsilon, betaGrowth, tstar, tref;
FILE * fpout = NULL;

/**
### Scriven solution */

double scriven_integrand (double x, void * params) {
  double b = *(double *)params;
  return exp (-sq(x) - 2.*epsilon*cube(b)/x)/sq(x);
}

double scriven_F (double s, double b) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_function F = {.function = scriven_integrand, .params = &b};
  double result, error;
  gsl_integration_qagiu (&F, s, 1e-12, 1e-8, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
  return result;
}

double scriven_phi (double b) {
  return 2.*cube(b)*exp (sq(b) + 2.*epsilon*sq(b))*scriven_F (b, b);
}

double scriven_beta (double Phi) {
  double a = 1e-3, b = 10.;
  for (int k = 0; k < 100; k++) {
    double c = 0.5*(a + b);
    if (scriven_phi (c) > Phi)
      b = c;
    else
      a = c;
  }
  return 0.5*(a + b);
}

double scriven_massfrac (double r, double time) {
  if (r < 2.*betaGrowth*sqrt (Dl*time))
    return wint;
  double s = r/(2.*sqrt (Dl*time));
  return wbulk - (wbulk - wint)*scriven_F (s, betaGrowth)
    /scriven_F (betaGrowth, betaGrowth);
}

int main (void) {

  /**
  The gas phase contains only the transferring species, while the liquid
  contains the dissolved gas and water. Species with the same name are matched
  between the two phases. */

  NGS = 1, NLS = 2;
  gas_species = (char *[]){"CO2"};
  liq_species = (char *[]){"CO2", "H2O"};

  /**
  Material properties from Table 1 of [Gennari et al., 2022](#gennari2022). The
  diffusivity is increased to $2\times10^{-5}$ m²/s to speed up the mass
  transfer. The gas is pure, so its diffusivity plays no role. */

  rho1 = 1000., rho2 = 1.8;
  mu1 = 1.05e-3, mu2 = 1.46e-5;
  Dmix1 = Dl, Dmix2 = Dl;

  f.sigma = 0.072;

  /**
  The molecular weights are equal, so that mole and mass fractions coincide,
  and the Henry constant becomes the interface mass fraction *YIntVal*. */

  MW1 = 1., MW2 = 1.;
  wint = rho2/(He*rho1);
  wbulk = ZETA*wint;
  YIntVal = wint;

  /**
  Growth constant and initial time. */

  epsilon = (rho1 - rho2)/rho1;
  double Phi = rho1/rho2*(wbulk - wint)/(1. - wint);
  betaGrowth = scriven_beta (Phi);
  tstar = sq(R0/(2.*betaGrowth))/Dl;
  tref = rho1*sq(2.*R0)/mu1;
  fprintf (stderr, "# beta = %g, t* = %g s, t*/tref = %g\n",
      betaGrowth, tstar, tstar/tref);

  /**
  Single velocity field, with the phase tracers advected consistently with
  the volume fraction. The system is isothermal. */

  nv = 1;
  pcm.consistent = true;
  pcm.isothermal = true;

  L0 = 16.*R0;
  TOLERANCE = 1e-6;

  /**
  We run three levels of refinement, i.e. 4, 8, and 16 cells per initial
  bubble radius. */

  for (maxlevel = 6; maxlevel <= 8; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

/**
The default settings of [evaporation.h](../src/evaporation.h) are for
droplets. For a bubble, the expansion term is shifted to the gas, so that the
liquid carries the Stefan flow, and the interface regression uses the liquid
density. This event runs after all the *init* events. */

event end_init (i = 0) {
  pcm.shifting = SHIFT_TO_GAS;
  pcm.byrhogas = false;
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {
  fraction (f, -circle (x, y, R0));

  /**
  The phase fields are tracers, i.e. multiplied by the volume fraction of
  their phase. */

  scalar YLCO2 = liq->YList[0], YLH2O = liq->YList[1];
  scalar YGCO2 = gas->YList[0];
  foreach() {
    double w = scriven_massfrac (sqrt (sq(x) + sq(y)), tstar);
    YLCO2[] = f[]*w;
    YLH2O[] = f[]*(1. - w);
    YGCO2[] = 1. - f[];
  }

  YLCO2[top] = dirichlet (wbulk);
  YLCO2[right] = dirichlet (wbulk);
  YLH2O[top] = dirichlet (1. - wbulk);
  YLH2O[right] = dirichlet (1. - wbulk);

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  fpout = fopen (name, "w");
}

#if TREE
event adapt (i++) {
  scalar YLCO2 = liq->YList[0];
  adapt_wavelet_leave_interface ({YLCO2, u.x, u.y}, {f},
      (double[]){1e-2*(wbulk - wint), 1e-4, 1e-4}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-processing

The effective radius is computed from the volume of the half bubble. The time
is shifted by $t^*$. */

double effective_radius() {
  double vol = 0.;
  foreach (reduction(+:vol))
    vol += (1. - f[])*dv();
  return cbrt (3.*vol);
}

event output_data (t += 0.01) {
  double time = t + tstar;
  fprintf (fpout, "%g %g %g\n", time, effective_radius(),
      2.*betaGrowth*sqrt (Dl*time));
}

/**
At the end of each simulation, we write the relative error of the radius,
and the concentration profile along the diagonal. */

event profiles (t = TMAX) {
  fclose (fpout);

  double R = effective_radius(), Rth = 2.*betaGrowth*sqrt (Dl*(t + tstar));
  double err = fabs (R - Rth)/Rth;
  fprintf (stderr, "%d %.4g\n", 1 << maxlevel, err);

  char name[80];
  sprintf (name, "Profile-%d", maxlevel);
  FILE * fp = fopen (name, "w");
  scalar YLCO2 = liq->YList[0];
  for (double r = 0.; r < 0.9*L0; r += L0/500.) {
    double xp = r*cos (pi/4.), yp = r*sin (pi/4.);
    double ff = interpolate (f, xp, yp);
    double w = interpolate (YLCO2, xp, yp);
    fprintf (fp, "%g %g %g\n", r, (ff > 1e-6) ? w/ff : wint,
        scriven_massfrac (r, t + tstar));
  }
  fclose (fp);
}

/**
## Results

~~~gnuplot Evolution of the bubble radius
reset
set xlabel "t [s]"
set ylabel "Bubble Radius [m]"
set key top left
set size square
set grid

plot "OutputData-8" every 20 u 1:3 w p t "Analytic", \
     "OutputData-6" u 1:2 w l lw 2 t "LEVEL 6", \
     "OutputData-7" u 1:2 w l lw 2 t "LEVEL 7", \
     "OutputData-8" u 1:2 w l lw 2 t "LEVEL 8"
~~~

The relative error of the radius at the final time, as a function of the
number of grid points $N$ per side of the domain.

~~~gnuplot Convergence of the bubble radius
reset
set xlabel "N"
set ylabel "Relative Error"

set xr[2**5:2**9]
set size square
set grid

set logscale x 2
set logscale y

f(x) = a*x**-b
fit f(x) "log" u 1:2 via a,b

ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", a, b)

plot "log" w p pt 8 title "Results", \
     f(x) w l lw 2 title ftitle(a, b)
~~~

~~~gnuplot Concentration profile
reset
set xlabel "Radius [m]"
set ylabel "Mass Fraction [-]"

set key bottom right
set size square
set grid

plot "Profile-8" every 10 u 1:3 w p t "Analytic", \
     "Profile-6" u 1:2 w l lw 2 t "LEVEL 6", \
     "Profile-7" u 1:2 w l lw 2 t "LEVEL 7", \
     "Profile-8" u 1:2 w l lw 2 t "LEVEL 8"
~~~

## References

~~~bib
@article{gennari2022,
  title={A phase-change model for diffusion-driven mass transfer problems in
    incompressible two-phase flows},
  author={Gennari, Gabriele and Jefferson-Loveday, Richard and Pickering,
    Stephen J},
  journal={Chemical Engineering Science},
  volume={259},
  pages={117791},
  year={2022},
  publisher={Elsevier}
}

@article{scriven1959,
  title={On the dynamics of phase growth},
  author={Scriven, LE},
  journal={Chemical Engineering Science},
  volume={10},
  number={1-2},
  pages={1--13},
  year={1959},
  publisher={Elsevier}
}
~~~
*/
