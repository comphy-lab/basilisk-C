/**
# Microgravity Droplet Evaporation

Evaporation of a droplet in microgravity conditions with variable material
properties.  This test case is sufficiently general to simulate different
initial conditions by overwriting the compilation macros defined in
[defaultvars.h](defaultvars.h). In this case, we simulate the same test case in
[c7pathak.c](c7pathak.c) which consists in the evaporation of pure n-heptane in
a nitrogen environment. By using the same input data we can test the effect of
variable material properties on the droplet consumption dynamics.

<div class="message"> Note that [burningdroplet.c](burningdroplet.c) provides
the same functionality, but with improved generality and flexibility. It should
be preferred.</div>
*/

/**
## Simulation Setup

We use the centered Navier--Stokes equations solver with volumetric source in
the projection step. The phase change is directly included using the evaporation
module, which sets the best (default) configuration for evaporation problems.
Many features of the phase change (evaporation) model can be modified directly
in this file without changing the source code, using the phase change model
object `pcm`. Compiling with `-DJUMP=1` changes the Navier--Stokes solver to the
velocity-jump formulation, which employs a GFM approach to set the interface
velocity jump. The density field is filtered in order to reduce convergence
issues related with strong density ratios. We include variable material properties
computed using OpenSMOKE++ or Cantera. */

#include "axi.h"
#if JUMP
# include "navier-stokes/velocity-jump.h"
#else
# include "navier-stokes/low-mach.h"
#endif
#define FILTERED 1
#define P_ERR 0.1
#include "two-phase-varprop.h"
#if USE_OPENSMOKE
# include "opensmoke/properties.h"
#elif USE_CANTERA
# include "cantera/properties.h"
#endif
#include "two-phase.h"
#include "tension.h"
//#include "recoil.h"
#include "evaporation.h"
#include "defaultvars.h"
#include "view.h"

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
### Variable properties functions

We declare variable properties functions for the liquid phase. This is not
necessary using OpenSMOKE++, but it is required by the Cantera interface, since
liquid phase laws are not implemented yet. */

double liqprop_density_heptane (void * p) {
  double a = 5.2745973, b = 0.07741, c = 557.342, d = 0.13673;
  ThermoState * ts = p;
  double Teff = min (ts->T, 0.999*c);
  return a / (pow (b, 1. + pow (1. - Teff / c, d)));
}

/**
### Simulation Data

We declare the maximum and minimum levels of refinement, the initial radius and
diameter, and the radius from the numerical simulation, and additional data for
post-processing. */

int maxlevel, minlevel = 2;
double D0 = DIAMETER, effective_radius0 = 1., d_over_d02 = 1., tad = 0.;
double volumecorr = 0., volume0 = 0., d_over_d02_stop = 1., mLiq0;
bool restored = false;

int main (void) {

  /**
  We set the kinetics folder, which defines the species of the simulation and 
  their properties. It is used by the kinetics library for the calculation
  of thermodynamic and transport properties. */

#if TWO_PHASE_VARPROP
  kinetics (TOSTRING(KINFOLDER), &NGS);
  kinetics_liquid (TOSTRING(KINFOLDER), &NLS);
  properties_liquid (TOSTRING(LIQFOLDER));
#else
  NGS = 2, NLS = 1;
#endif

  /**
  We set additional simulation properties. Be careful, these are "dummy"
  properties which are overwritten by the variable-properties formulation
  already at the first iteration. To start the simulation without any
  problem of division by zero, the viscosity value should be non-null. If
  a function in [ThermoProps](../src/variable-properties.h) is NULL, the
  value initialized here is used. */

  rho1 = 626.7; rho2 = 17.51;
  mu1 = 1.e-3; mu2 = 1.e-5;
  Dmix1 = 0., Dmix2 = 6.77e-7;
  lambda1 = 0.1121, lambda2 = 0.04428;
  cp1 = 2505., cp2 = 1053.;
  dhev = 3.23e5;

  /**
  We set the thermodynamic pressure in SI units (Pa), and the initial
  temperatures of the system. */

  Pref = PRESSURE*101325.;
  TG0 = TEMPERATURE;
  TL0 = TEMPERATURE_DROPLET;

  /**
  We solve two different sets of Navier--Stokes equations according with the
  double pressure velocity coupling approach, or with the velocity jump
  formulations.  In either case we need two different velocity fields. Here
  we can overwrite some phase change model `pcm` property. The emissivity
  controls the interface radiation. */

  nv = 2;
  pcm.emissivity = EMISSIVITY;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  double RR = ENVIRONMENT;
  L0 = 0.5*RR*D0;

  /**
  We change the surface tension coefficient. */

  f.sigma = SIGMA;

  /**
  We run the simulation at different maximum
  levels of refinement. */

  maxlevel = MAXLEVEL;
  init_grid (1 << min (maxlevel, 9));
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

event init (i = 0) {
  if (!restore (file = "restart")) {
#if TREE
    refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);
#endif
    fraction (f, circle (x, y, 0.5*D0));

    /**
    We compute initial variables useful for post-processing. */

    effective_radius0 = pow (3.*statsf(f).sum, 1./3.);
    volume0 = 4./3.*pi*pow (effective_radius0, 3.);

    scalar rhol = liq->rho;
    foreach (reduction(+:mLiq0))
      mLiq0 += rhol[]*f[]*dv();
  }
  else
    restored = true;

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  /**
  We set the properties of the system. */

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = (double[]){1.};
  tsg.T = TG0, tsg.P = Pref, tsg.x = (double[]){0.,1.};

  phase_set_thermo_state (liq, &tsl, force = !restored);
  phase_set_thermo_state (gas, &tsg, force = !restored);

  phase_set_properties (liq, MWs = (double[]){100.2});
  phase_set_properties (gas, MWs = (double[]){100.2,29.});

#if OPENSMOKE
  antoine = &opensmoke_antoine;
#else
  tp1.rhov = liqprop_density_heptane;

  scalar lfuel = liq->YList[0];
  lfuel.antoine = antoine_heptane;
#endif

  /**
  On the top and right sides we set Dirichlet boundary conditions for the
  temperature and mass fraction fields. */

  scalar fuel  = gas->YList[0];
  scalar inert = gas->YList[1];
  scalar TG    = gas->T;

  fuel[top] = dirichlet (0.);
  fuel[right] = dirichlet (0.);

  inert[top] = dirichlet (1.);
  inert[right] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({Y,T,u.x,u.y}, {f},
      (double[]){YTOL,TTOL,UTOL,UTOL}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Logger

We output the dimensionless droplet radius. */

event logger (t += 0.000016) {
  double R = pow (3.*statsf(f).sum, 1./3.);
  fprintf (stderr, "%d %g %.5g\n", i, t, R / (0.5*D0));
}

/**
### Output Files

We write on a file the squared diameter decay and the
dimensionless time. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double effective_radius = pow (3.*statsf(f).sum, 1./3.);
  double d_over_d02_old = d_over_d02;
  double tad_old = tad;

  d_over_d02 = sq (effective_radius / effective_radius0);
  d_over_d02_stop = sq (effective_radius / (0.5*D0));
  tad = t/sq(D0*1e3);

  /**
  The vaporization rate is computed according to the formula
  in Liu & Avedisian, 2011, pag. 777 bottom. */

  double kv = 0.;
  if (i > 1)
    kv = fabs ((d_over_d02 - d_over_d02_old)/(tad - tad_old));

  double mLiq = 0.;
  scalar rhol = liq->rho;
  foreach(reduction(+:mLiq))
    mLiq += rhol[]*f[]*dv();

  /**
  We compute and print additional useful average quantities. */

  double TIntAvg = avg_interface (gas_int->T, f, tol=0.1);
  double YIntAvg = avg_interface (gas_int->YList[0], f, tol=0.1);

  int counter = 0;
  double TDropAvg = 0.;
  scalar TL = liq->T;
  foreach(reduction(+:TDropAvg) reduction(+:counter)) {
    if (f[] > 1.-F_ERR) {
      counter++;
      TDropAvg += TL[];
    }
  }
  TDropAvg = (counter > 0.) ? TDropAvg/counter : 0.;

  fprintf (fp, "%g %g %g %g %g %g %g %g %g %g\n", t, tad, effective_radius0,
      effective_radius, d_over_d02, mLiq, kv, TIntAvg, YIntAvg, TDropAvg);
}

/**
### Movie

We write the animation with the evolution of the
n-heptane mass fraction, the interface position
and the temperature field. */

#if MOVIE
event movie (t += MOVIE_EVERY) {
  clear();
  box();
  view (fov = 3, samples = 2);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = TL0, max = statsf(T).max, linear = true);
  save ("temperature.mp4");

  clear();
  box();
  view (fov = 3, samples = 2);
  draw_vof ("f", lw = 1.5);
  squares ("Y", min = 0., max = 1., linear = true);
  save ("fuel.mp4");
}
#endif

/**
### Snapshots

Output dump files for restore or post-processing. */

#if DUMP
event snapshots (t += DUMP_EVERY) {
  char name[80];
  sprintf (name, "snapshots-%g", t);
  dump (name);
}
#endif

/**
### Stopping Condition

We stop the simulation when the droplet is almost fully consumed. */

event stop (i++) {
  if (d_over_d02_stop <= MAX_DD02)
    return 1;
}

/**
We run the simulation for long time, which is not reached because
the stopping condition on the droplet diameter is reached first. */

#if BASILISK_SERVER
event end (t = 1.6e-4);
#else
event end (t = 50.);
#endif

/**
## Results

We compare the constant and variable properties results. With variable
properties we can capture phenomena such as the initial thermal expansion, and
the different steady vaporization rate constant. Using OpenSMOKE++ the gap
between the two models is larger because it implements variable liquid
properties. With Cantera we should introduce more functions for the liquid
properties.

~~~gnuplot Square diameter decay
set xlabel "time [s]"
set ylabel "(D/D_0)^2 [-]"
set grid

p "OutputData-6" u 1:($4/2.5e-6)**2 w l t "Variable properties", \
  "../c7pathak/OutputData-6" u 1:($3/2.5e-6)**2 w l t "Constant properties"
~~~

~~~gnuplot
set xlabel "time [s]"
set ylabel "Interface temperature [K]"
set grid
set yr[380:]

p "OutputData-6" u 1:8 w l t "Variable properties", \
  "../c7pathak/OutputData-6" u 1:5 w l t "Constant properties"
~~~
*/
