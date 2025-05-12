/**
# Thermal Expansion of a Liquid Droplet

In this test case, we simulate the thermal expansion of an initially cold
n-heptane droplet in a hot isothermal environment, neglecting the phase
change. The aim of this simulation is to test the convergence of the
volumetric source term due to density changes, using the
[multicomponent.h](/sandbox/ecipriano/src/multicomponent-varprop.h) solver with variable
properties.

![Evolution of the temperature field](expansion/movie.mp4)
*/

/**
## Phase Change Setup

We define the number of gas and liquid species, and we set
the initial composition of the two phases, which does not change in time. 
By setting the thermodynamic equilibrium `Keq` to 0 we suppress the phase
change. */

#define NGS 2
#define NLS 1

char* gas_species[NGS] = {"NC7H16", "N2"};
char* liq_species[NLS] = {"NC7H16"};
char* inert_species[1] = {"N2"};
double gas_start[NGS] = {0., 1.};
double liq_start[NLS] = {1.};
double inKeq[NLS] = {0.};

/**
We declare all the properties necessary for the multicomponent phase
change model. The properties are set to the initial bulk values and they are
overwritten by the variable-properties formulation, which computes all the
physical properties as a function of the thermodynamic state of the mixture. */

double inDmix1[NLS] = {0.};
double inDmix2[NGS] = {9.19211e-07, 3.28466e-06};
double lambda1 = 0.124069;
double lambda2 = 0.0295641;
double dhev = 364482;
double cp1 = 2244.92;
double cp2 = 1041.52;

/**
We set the initial temperature of the liquid and of the gas phase. The gas
phase temperature can be changed by the user through a compilation variable,
for example by setting `-DTEMPERATURE=400` */

#ifndef TEMPERATURE
# define TEMPERATURE 350
#endif
double TL0 = 300.;
double TG0 = TEMPERATURE;

/**
We solve the temperature field, and we reduce the tolerance for the
calculation of the variable properties. The interfacial
temperature is not computed from the jump conditon, it is just set to the gas
phase temperature value, in order to keep the gas phase temperature constant.
The advection term of the momentum equation is solved in non-conservative form 
by setting the variable `NO_ADVECTION_DIV` to 1. */

#define SOLVE_TEMPERATURE
#define NO_ADVECTION_DIV 1
#define T_PROP 1.e-6
#define FIXED_INTERFACE_TEMPERATURE TG0

/**
## Simulation Setup

In this simulation, the choice of the Navier-Stokes equations solver is not
important, because there is no velocity jump due to the phase change.
We can choose between the centered solver with expansion source term, and
the centered solver with velocity jump.
The [two-phase.h](/src/two-phase.h) solver is extended to variable physical
properties by including a policy for the calculation of such properties. In
this case, we use correlations implemented in OpenSMOKE++ (refer to
[Cipriano et al., 2024](#cipriano2024) for details).
The solution of the temperature field is performed by the multicomponent
phase change model.
*/

#include "axi.h"
#if JUMP
# include "navier-stokes/velocity-jump.h"
#else
# include "navier-stokes/centered-evaporation.h"
# include "navier-stokes/centered-doubled.h"
#endif
#if BASILISK_SERVER
# include "basilisk-properties.h"
#else
# include "opensmoke-properties.h"
#endif
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "multicomponent-varprop.h"
#include "view.h"

/**
### Boundary conditions

We initialize the droplet at the lower left corner of the domain, and we
exploit axial symmetry. We set symmetry boundary conditions on the sides
in contact with the droplet (left, bottom), and outflow boundary
conditions on the other sides of the domain. */

#if JUMP
u1.n[top] = neumann (0.);
u1.t[top] = neumann (0.);
u2.n[top] = neumann (0.);
u2.t[top] = neumann (0.);
p[top] = dirichlet (0.);
ps[top] = dirichlet (0.);
pg[top] = dirichlet (0.);

u1.n[right] = neumann (0.);
u1.t[right] = neumann (0.);
u2.n[right] = neumann (0.);
u2.t[right] = neumann (0.);
p[right] = dirichlet (0.);
ps[right] = dirichlet (0.);
pg[right] = dirichlet (0.);
#else
u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
uext.n[top] = neumann (0.);
uext.t[top] = neumann (0.);
pext[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);
#endif

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial diameter (1 mm), and additional data for
post-processing. */

int maxlevel, minlevel = 2;
double D0 = 1e-3, effective_radius0, d_over_d02 = 1.;
double volume0 = 0., rho0 = 0., rhof = 0.;

int main (void) {
  
  /**
  We set the kinetics folder, which defines the species
  of the simulation, and it is used by OpenSMOKE++ for the
  calculation of the thermodynamic and transport properties.
  The path is relative to `OpenSMOKEppInterface/kinetics`. */

#if !BASILISK_SERVER
  kinfolder = "evaporation/n-heptane-in-nitrogen";
#endif

  /**
  We set additional simulation properties. Be careful, these are "dummy"
  properties which are overwritten by the variable-properties formulation
  already at the first iteration. To start the simulation without any
  problem of division by zero, the viscosity value should be non-null.
  All the other properties can be null if proper thermodynamic rules are
  provided for all those properties. */

  rho1 = 681.042; rho2 = 9.75415;
  mu1 = 0.00037446; mu2 = 2.02391e-05;
  
  /**
  The thermodynamic pressure is set to 10 atm, we consider a `Pref` value
  in SI units (Pa). */
  
  Pref = 10*101325.;

  /**
  We change the dimensions of the domain as a function of the initial
  diameter of the droplet. */

  L0 = 1.5*D0;

  /**
  We use a constant surface tension value. */

  f.sigma = 0.03;

  /**
  We run the simulation at different maximum levels of refinement. */

  for (maxlevel = 5; maxlevel <= 5; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field according to `D0`. */

event init (i = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  /**
  We compute initial variables useful for post-processing. */

  effective_radius0 = pow (3.*statsf(f).sum, 1./3.);
  volume0 = 4./3.*pi*pow (effective_radius0, 3.);
  
  foreach (reduction(+:mLiq0))
    mLiq0 += rho1v[]*f[]*dv();

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

#if BASILISK_SERVER
  inMW[0] = 100.; inMW[1] = 29.;

  /**
  We set the functions for variable liquid density. Note
  that if using OpenSMOKE++, every variable changes as
  explained in [Cipriano et. al, 2024](#cipriano2024b). */

  tp1.rhov = liqprop_density_heptane;
#else
  for (int jj=0; jj<NGS; jj++)
    inMW[jj] = OpenSMOKE_MW (jj);
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

#if PRINT_PROPERTIES
  FILE * fp1 = fopen ("liqprop", "w");
  FILE * fp2 = fopen ("gasprop", "w");

  ThermoState ts1;
  ts1.T = TL0;
  ts1.P = Pref;
  ts1.x = liq_start;

  ThermoState ts2;
  ts2.T = TG0;
  ts2.P = Pref;
  ts2.x = gas_start;

  print_thermoprop (&tp1, &ts1, NLS, fp=fp1);
  print_thermoprop (&tp2, &ts2, NGS, fp=fp2);

  fclose (fp1);
  fclose (fp2);
#endif

  /**
  On the top and right sides we set Dirichlet boundary conditions
  for the temperature and mass fraction fields. */

  scalar NC7H16 = YGList[0];
  scalar N2    = YGList[1];

  NC7H16[top] = dirichlet (0.);
  NC7H16[right] = dirichlet (0.);

  N2[top] = dirichlet (1.);
  N2[right] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

/**
We adapt the grid according to the temperature and the velocity fields. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({T,u.x,u.y}, {f},
      (double[]){1.e-2,1.e-1,1.e-1}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Output Files

We write on a file squared diameter, the exact droplet diameter at
steady-state, and the ratio between the current and the initial liquid
mass values. */

event output_data (i += 50) {
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
      d_over_d02, mLiq/mLiq0, relerr, exact_d_over_d02);
}

/**
### Movie

We write images with the temperature and the divergence of the liquid
velocity at time 0.1 s. */

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
    drhodtplot[] = drhodtext[]*f[];

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
}

/**
We write a video with the evolution of the temperature field. */

event movie (t += 0.02; t <= 3) {
  clear();
  box();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f");
  squares ("T", min = TL0, max = statsf(T).max, linear = true);
  save ("movie.mp4");
}

/**
## Results

Trend of the square droplet diameter in time. The droplets expands differently
depending on the gradient between the interface and the gas phase temperature.
We just run the lowest level of refinement for fast simulations on the Basilisk
server. The simulation converges to the steady state values with increasing
mesh resolution.

~~~gnuplot Expansion of the Square Diameter
set border linewidth 2
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

plot basedir350."OutputData-5" u 2:7 every 100 w p ps 0.8 lc 1 pt 6 dt 2 t "Steady State", \
     basedir350."OutputData-5" u 2:4 w l lw 2 lc 1 dt 4 t "LEVEL 5", \
     basedir375."OutputData-5" u 2:7 every 100 w p ps 0.8 lc 2 pt 6 dt 2 notitle, \
     basedir375."OutputData-5" u 2:4 w l lw 2 lc 2 dt 4 notitle, \
     basedir400."OutputData-5" u 2:7 every 100 w p ps 0.8 lc 3 pt 6 dt 2 notitle, \
     basedir400."OutputData-5" u 2:4 w l lw 2 lc 3 dt 4 notitle
~~~
*/
