/**
# Expansion of a Liquid Droplet

In this test case, we study the thermal expansion of a liquid droplet of
n-heptane, at different ambient temperatures and at different levels of
refinement. The aim is to evaluate the convergence of the variable material
properties solver. Evaporation is neglected in order to focus exclusively on the
thermal expansion effects.

![Evolution of the temperature field](expansion/movie.mp4)
*/

/**
## Simulation Setup

For variable material properties we need a Navier--Stokes solver with volumetric
sources, whcih account for the expansion (or compression) due to thermal or
composition effects.  The [two-phase-varprop.h](../src/two-phase-varprop.h)
overwrite the one-field density and viscosity in order to include the effect of
variable phase properties. The OpenSMOKE++ library is used to calculate the
phase properties, while the phase change model is used to resolve the chemical
species and temperature fields and for the interface boundary conditions.
*/

#include "axi.h"
#if JUMP
# include "navier-stokes/velocity-jump.h"
#else
# include "navier-stokes/low-mach.h"
#endif
#define P_ERR 1.e-6
#include "two-phase-varprop.h"
#include "basilisk-properties.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "balances/two-phase.h"
#include "view.h"

/**
The following data can be overwritten during compilation. */

#ifndef TEMPERATURE
# define TEMPERATURE 350
#endif
#ifndef PRESSURE
# define PRESSURE 10
#endif

/**
We define the function/laws for updating material properties. In this case, a
function for the calculation of liquid phase density. */

#if BASILISK_PROPERTIES
double liqprop_density_heptane (void * p) {
  double a = 5.2745973, b = 0.07741, c = 557.342, d = 0.13673;
  ThermoState * ts = p;
  double Teff = min (ts->T, 0.999*c);
  return a / (pow (b, 1. + pow (1. - Teff / c, d)));
}
#endif

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right sides of the domain. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

/**
### Simulation Data

We declare the maximum and minimum levels of refinement, the initial radius and
diameter, and the radius from the numerical simulation, and additional data for
post-processing. */

int maxlevel, minlevel = 2;
double D0 = 1e-3, effective_radius0, d_over_d02 = 1., mLiq0;
double volumecorr = 0., volume0 = 0., rho0 = 0., rhof = 0.;

int main (void) {
  /**
  We set the kinetics folder, which defines the species of the simulation, and
  it is used by OpenSMOKE++ for the calculation of the thermodynamic and
  transport properties. */

#if defined(TWO_PHASE_VARPROP) && !defined(BASILISK_PROPERTIES)
  kinetics ("evaporation/n-heptane-in-nitrogen", &NGS);
  kinetics_liquid ("evaporation/n-heptane-in-nitrogen", &NLS);
  properties_liquid ("LiquidProperties");
#else
  NGS = 2, NLS = 1;
#endif

  /**
  We set the material properties for this simulation. They are used just in case
  OpenSMOKE++ is not installed in our system (i.e. on the basilisk server) and
  we do not want/need to change every single property.  Alternatively, if we use
  OpenSMOKE++ (or another library) we don't need it. */

  rho1 = 681.042; rho2 = 9.75415;
  mu1 = 0.00037446; mu2 = 2.02391e-05;
  Dmix1 = 0., Dmix2 = 9.19211e-07;
  lambda1 = 0.124069, lambda2 = 0.0295641;
  cp1 = 2244.92, cp2 = 1041.52;
  dhev = 364482;

  Pref = PRESSURE*101325.;
  TG0 = TEMPERATURE, TL0 = 300, TIntVal = TG0;

  /**
  The interface is considered isothermal and its temperature is `TIntVal`. */

  nv = 2;
  pcm.isothermal_interface = true;
  pcm.divergence = true;

  /**
  We change the dimension of the domain as a function of the initial diameter of
  the droplet. */

  L0 = 1.5*D0;

  /**
  We change the surface tension coefficient. The value is reduced to a very small
  surface tension compared to the value reported in the paper, because we want to
  reduce the simulation time on the Basilisk server for quick testing. */

  //f.sigma = 0.03;
  f.sigma = 0.0001;

  /**
  We run the simulation at different maximum levels of refinement. */

  for (maxlevel = 5; maxlevel <= 7; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the initial radius of the
droplet. We don't use the value D0 because for small errors of initialization
the squared diameter decay would not start from 1. */

event init (i = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  /**
  We compute initial variables useful for post-processing. */

  effective_radius0 = pow (3.*statsf(f).sum, 1./3.);
  volume0 = 4./3.*pi*pow (effective_radius0, 3.);

  scalar rhol = liq->rho;
  foreach (reduction(+:mLiq0))
    mLiq0 += rhol[]*f[]*dv();

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = (double[]){1.};
  tsg.T = TG0, tsg.P = Pref, tsg.x = (double[]){0.,1.};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  phase_set_properties (liq, MWs = (double[]){100.2});
  phase_set_properties (gas, MWs = (double[]){100.2,29.});

#if BASILISK_PROPERTIES
  tp1.rhov = liqprop_density_heptane;
#endif

  /**
  We compute variables for the analytical steady-state solution. */

  ThermoState ts0;
  ts0.T = TL0;
  ts0.P = Pref;
  ts0.x = (double[]){1.};

  ThermoState tsf;
  tsf.T = TG0;
  tsf.P = Pref;
  tsf.x = (double[]){1.};

  rho0 = tp1.rhov (&ts0);
  rhof = tp1.rhov (&tsf);

  /**
  We add boundary conditions for the scalar fields. In particular, we fix the
  temperature and composition on the top and right sides of the domain. */

  scalar NC7H16 = gas->YList[0], N2 = gas->YList[1];
  scalar TG = gas->T;

  NC7H16[top] = dirichlet (0.);
  NC7H16[right] = dirichlet (0.);

  N2[top] = dirichlet (1.);
  N2[right] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

/**
We adapt the grid according to the mass fractions of the mass fraction of
n-heptane, the temperature, and the velocity field. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({T,u.x,u.y}, {f},
      (double[]){1e0,1.e-1,1.e-1}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Output Files

We write on a file the squared diameter decay and the dimensionless time. */

event output_data (t += 0.02) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double effective_radius = pow (3.*statsf(f).sum, 1./3.);
  d_over_d02 = sq (effective_radius / effective_radius0);

  double mLiq = 0.;
  foreach(reduction(+:mLiq))
    mLiq += rho1v[]*f[]*dv();

  double exact_volume = volume0*rho0/rhof;
  double exact_radius = pow (3./4./pi*exact_volume, 1./3.);
  double exact_d_over_d02 = sq (exact_radius/effective_radius0);
  double relerr = fabs (effective_radius - exact_radius)/effective_radius;

  fprintf (fp, "%g %g %g %g %g %g %g\n", t, t/sq(D0*1e3), effective_radius,
      d_over_d02, mLiq/mLiq0, relerr, exact_d_over_d02), fflush (fp);
}

/**
### Logger

We output the total liquid volume in time (for testing). */

event logger (t += 0.2) {
  fprintf (stderr, "%d %.2f %.3g\n", i, t, statsf (f).sum);
}

/**
### Movie

We output pictures of the temperature and velocity divergence inside the
droplet. */

event pictures (t = 0.1) {
  clear();
  view (ty = -0.5);
  draw_vof ("f", filled = -1, fc = {1.,1.,1.});
  squares ("T", min = TL0, max = TG0, linear = true);
  isoline ("T", n = 12., min = TL0, max = TG0);
  mirror ({1.,0.}) {
    draw_vof ("f", filled = -1, fc = {1.,1.,1.});
    squares ("T", min = TL0, max = TG0, linear = true);
    isoline ("T", n = 12., min = TL0, max = TG0);
  }
  save ("temperature.png");

  scalar drhodtplot[];
  foreach()
    drhodtplot[] = drhodt[]*f[];

  clear();
  view (ty = -0.5);
  draw_vof ("f", filled = -1, fc = {1.,1.,1.});
  squares ("drhodtplot", spread = -1, linear = true);
  isoline ("drhodtplot", n = 12.,
      min = statsf(drhodtplot).min, max = statsf(drhodtplot).max);
  mirror ({1.,0.}) {
    draw_vof ("f", filled = -1, fc = {1.,1.,1.});
    squares ("drhodtplot", spread = -1, linear = true);
    isoline ("drhodtplot", n = 12.,
        min = statsf(drhodtplot).min, max = statsf(drhodtplot).max);
  }
  save ("divergence.png");

  scalar un[];
  foreach()
    un[] = norm (u)*f[];
  clear();
  view (ty = -0.5);
  draw_vof ("f", filled = -1, fc = {1.,1.,1.});
  squares ("un", spread = -1, linear = true);
  isoline ("un", n = 12., min = statsf(un).min, max = statsf(un).max);
  mirror ({1.,0.}) {
    draw_vof ("f", filled = -1, fc = {1.,1.,1.});
    squares ("un", spread = -1, linear = true);
    isoline ("un", n = 12., min = statsf(un).min, max = statsf(un).max);
  }
  save ("velocity.png");
}

/**
### Movie

We write the animation with the evolution of the n-heptane mass fraction, the
interface position and the temperature field. */

event movie (t += 0.05; t <= 3) {
  if (TG0 == 350 && maxlevel == 7) {
    clear();
    box();
    view (tx = -0.5, ty = -0.5);
    draw_vof ("f");
    squares ("T", min = TL0, max = statsf(T).max, linear = true);
    save ("movie.mp4");
  }
}

/**
## Results

Trend of the square droplet diameter in time. The droplets expands differently
depending on the gradient between the interface and the gas phase temperature.
We just run the lowest level of refinement for fast simulations on the Basilisk
server. The simulation converges to the steady state values with increasing
mesh resolution.

~~~gnuplot Expansion of the Square Diameter
set grid
set key bottom right
set xlabel "t/D_0^2 [s/mm^2]"
set ylabel "(D/D_0)^2 [-]"
set xrange[0:3]
set yrange[1:1.12]
set size square

basedir350 = "../expansion/"
basedir375 = "../expansion-T375/"
basedir400 = "../expansion-T400/"

set label "∆T = 50 K"  at 2.45,1.05  left font ",11" tc rgb "black"
set label "∆T = 75 K"  at 2.45,1.078 left font ",11" tc rgb "black"
set label "∆T = 100 K" at 2.45,1.11  left font ",11" tc rgb "black"

plot basedir350."OutputData-5" u 2:7 every 10 w p ps 0.8 lc 1 pt 6 dt 2 t "Steady State", \
     basedir350."OutputData-5" u 2:4 w l lc 1 dt 4 t "LEVEL 5", \
     basedir350."OutputData-6" u 2:4 w l lc 1 dt 3 t "LEVEL 6", \
     basedir350."OutputData-7" u 2:4 w l lc 1 dt 1 t "LEVEL 7", \
     basedir375."OutputData-5" u 2:7 every 10 w p ps 0.8 lc 2 pt 6 dt 2 notitle, \
     basedir375."OutputData-5" u 2:4 w l lc 2 dt 4 notitle, \
     basedir375."OutputData-6" u 2:4 w l lc 2 dt 3 notitle, \
     basedir375."OutputData-7" u 2:4 w l lc 2 dt 1 notitle, \
     basedir400."OutputData-5" u 2:7 every 10 w p ps 0.8 lc 3 pt 6 dt 2 notitle, \
     basedir400."OutputData-5" u 2:4 w l lc 3 dt 4 notitle, \
     basedir400."OutputData-6" u 2:4 w l lc 3 dt 3 notitle, \
     basedir400."OutputData-7" u 2:4 w l lc 3 dt 1 notitle
~~~

~~~gnuplot Convergence rate
reset

stats "<tail -n 1 OutputData-5" u 6 nooutput name "T1_LEVEL5"
stats "<tail -n 1 OutputData-6" u 6 nooutput name "T1_LEVEL6"
stats "<tail -n 1 OutputData-7" u 6 nooutput name "T1_LEVEL7"

stats "<tail -n 1 ../expansion-T375/OutputData-5" u 6 nooutput name "T2_LEVEL5"
stats "<tail -n 1 ../expansion-T375/OutputData-6" u 6 nooutput name "T2_LEVEL6"
stats "<tail -n 1 ../expansion-T375/OutputData-7" u 6 nooutput name "T2_LEVEL7"

stats "<tail -n 1 ../expansion-T400/OutputData-5" u 6 nooutput name "T3_LEVEL5"
stats "<tail -n 1 ../expansion-T400/OutputData-6" u 6 nooutput name "T3_LEVEL6"
stats "<tail -n 1 ../expansion-T400/OutputData-7" u 6 nooutput name "T3_LEVEL7"

set print "errors"

print sprintf ("350K %d %.12f", 2**5, T1_LEVEL5_mean)
print sprintf ("350K %d %.12f", 2**6, T1_LEVEL6_mean)
print sprintf ("350K %d %.12f", 2**7, T1_LEVEL7_mean)

print sprintf ("375K %d %.12f", 2**5, T2_LEVEL5_mean)
print sprintf ("375K %d %.12f", 2**6, T2_LEVEL6_mean)
print sprintf ("375K %d %.12f", 2**7, T2_LEVEL7_mean)

print sprintf ("400K %d %.12f", 2**5, T3_LEVEL5_mean)
print sprintf ("400K %d %.12f", 2**6, T3_LEVEL6_mean)
print sprintf ("400K %d %.12f", 2**7, T3_LEVEL7_mean)

unset print

reset
set xlabel "N"
set ylabel "Relative Error"

set xr[2**4:2**8]
set size square
set grid

set logscale x 2
set logscale y

f(x) = a*x**-b

# fit for 350K
a = 1; b = 1
fit f(x) "<grep 350K errors" u 2:3 via a,b
a350 = a; b350 = b

# fit for 375K
a = 1; b = 1
fit f(x) "<grep 375K errors" u 2:3 via a,b
a375 = a; b375 = b

# fit for 400K
a = 1; b = 1
fit f(x) "<grep 400K errors" u 2:3 via a,b
a400 = a; b400 = b

ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)

plot "<grep 350K errors" u 2:3 w p pt 4 lc 1 title "350K data", \
     "<grep 375K errors" u 2:3 w p pt 6 lc 2 title "375K data", \
     "<grep 400K errors" u 2:3 w p pt 8 lc 3 title "400K data", \
     a350*x**-b350 w l lc 1 title sprintf("350K: %s", ftitle(a350,b350)), \
     a375*x**-b375 w l lc 2 title sprintf("375K: %s", ftitle(a375,b375)), \
     a400*x**-b400 w l lc 3 title sprintf("400K: %s", ftitle(a400,b400))
~~~
*/

