/**
# Evaporation of a Pure Droplet in Buoyancy-Driven Flows

This module reports the simulation setup for a generic pure droplet
at different gravity conditions. This configuration has been used in
several experimental works, by suspending the droplet in a static
configuration by means of a solid fiber. This module mimick the
experimental setup by combining a droplet suspension strategy with
the interface-resolved phase change model and the variable properties
formulation, which allows buoyancy-driven flows to be resolved directly.
These phenomena promote non-radial gas-phase convective fluxes, and
liquid phase internal recirculation, which break the spherical symmetry
of the evaporation process in microgravity conditions.
The simulation setup can be used with any gravity value, not only in
normal gravity. This feature allows small residual gravity values to be
simulated, making this simulation setup useful to drive the experimental
investigation also for droplets in microgravity.

This simulation includes variable thermodynamic and transport properties,
and the interface radiation contribution. We neglect the heat transfer
from the suspending solid fiber, which may affect experimental data if
the fiber is large. The simulation is performed by simulating a droplet
at the edge of a square domain, exploting the axial-symmetry, which is
acceptable for small Reynolds number.

<div class="message">
Note that [burningdroplet.c](burningdroplet.c) provides the same functionality,
but with improved generality and flexibility. It should be preferred.</div>
*/

/**
## Default Simulation Data

The simulation setup is independent from the liquid under investigation,
and we may want to perform the same simulation studying the sensitivity
of the numerical results to the operative conditions (temperature,
pressure, ecc...) and to the liquid fuel or the ambient conditions.
To avoid code duplication, we set a bunch of default compiler variables,
which are overwritten in the `Makefile` to create different cases with
different operative conditions. The default properties are collected in the
file [defaultvars.h](defaultvars.h).

By default, in this test case we consider the properties of a n-heptane liquid,
evaporating in pure nitrogen environment. We considers a 1 mm droplet, without
interface radiation, in a 773 K environment. The simulation runs until 0.5 s.
During this transient, we observe the heating of the liquid phase and the
formation of a downward wake (P = 10 atm). The simulation is stopped before the
complete consumption of the droplet in order to reduce the simulation time on
the basilisk server as much as possible.

![Temperature field](normalgravity/temperature.mp4)
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
computed using OpenSMOKE++ or Cantera. Using
[gravity.h](/sandbox/ecipriano/src/gravity.h) simplifies the boundary condition
for pressure, analogously to the [reduced.h](/src/reduced.h) approach but
considering variable liquid and gas phase densities.
*/

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
#include "pinning.h"
#include "two-phase.h"
#include "tension.h"
#include "gravity.h"
//#include "recoil.h"
#include "evaporation.h"
#include "defaultvars.h"
#include "view.h"

/**
### Boundary conditions

We initialize the droplet at the lower edge of the domain, and we exploit axial
symmetry. We set no-slip boundary conditions on the side in contact with the
droplet (left), and outflow boundary conditions on the other sides of the
domain. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

u.n[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);
p[bottom] = neumann (0.);
uf.n[bottom] = 0.;
uf.t[bottom] = 0.;

/**
### Simulation Data

We declare the maximum and minimum levels of refinement, the initial droplet
diameter, and additional data for post-processing. We also transport a scalar
tracer `tr` which can be useful to visualize liquid internal recirculation. */

int maxlevel, minlevel = 2;
double D0 = DIAMETER, effective_radius0, mLiq0;
double effective_radius = 0.5*DIAMETER, d_over_d02 = 1., tad = 0.;
double volumecorr = 0., trmin = 0., trmax = 0.;
bool restored = false;

scalar tr[];

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

  rho1 = 681.042; rho2 = 9.75415;
  mu1 = 0.00037446; mu2 = 2.02391e-05;
  Dmix1 = 0., Dmix2 = 1e-4;
  lambda1 = 0.124069, lambda2 = 0.0295641;
  cp1 = 2244.92, cp2 = 1041.52;
  dhev = 364482;

  /**
  We set the thermodynamic pressure in SI units (Pa), and the initial
  temperatures of the system. */

  Pref = PRESSURE*101325.;
  TG0 = TEMPERATURE, TL0 = TEMPERATURE_DROPLET;

  /**
  We solve two different sets of Navier--Stokes equations according with the
  double pressure velocity coupling approach, or with the velocity jump
  formulations.  In either case we need two different velocity fields. Here
  we can overwrite some phase change model `pcm` property. The emissivity
  controls the interface radiation. */

  nv = 2;
  pcm.emissivity = EMISSIVITY;

  /**
  We change the dimensions of the domain as a function of the initial
  diameter of the droplet. The domain is large with respect to the
  droplet diameter, in order to avoid the influence of the domain
  boundaries on the droplet evaporation dynamics. */

  double RR = ENVIRONMENT;
  L0 = 0.5*RR*D0;

  /**
  We set the gravity contribution, and we shift the origin along the `y`
  coordinate, according to the radius of the solid fiber. */

  G.x = GRAVITY;
  double df = FIBER*D0;
  X0 = -0.5*L0, Y0 = 0.5*df;

  /**
  We use a constant surface tension value, and we transport `tr` as
  a VOF tracer. */

  f.sigma = 0.03;
  f.tracers = {tr};

  /**
  We run the simulation at different maximum levels of refinement. */

  maxlevel = MAXLEVEL;

  /**
  The pinning point `ap` is adjusted according to the value of `Y0`.
  The pinning center `ac` is also adjusted depending on the grid. */

  pinning.ap = sqrt (sq (0.5*D0) - sq (Y0));
  pinning.ac = pinning.ap - 2.*L0/(1 << maxlevel);

  init_grid (1 << min (maxlevel, 9));
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field according to `D0`,
after refining the region around the droplet to the maximum
level of refinement. */

event init (i = 0) {
  if (!restore (file = "restart")) {
    refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);
    fraction (f, circle (x, y, 0.5*D0));

    /**
    We compute initial variables useful for post-processing. */

    //volumecorr = 2.*pi*statsf(f).sum - (4./3.*pi*pow (0.5*D0, 3.));
    volumecorr = 0.;
    effective_radius0 = pow(3./4./pi*(2.*pi*statsf(f).sum - volumecorr), 1./3.);
    effective_radius = effective_radius0;

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
  On the edge in contact with the liquid droplet we set Neumann
  boundary conditions for the scalar fields, while we set Dirichlet
  boundary conditions for the other sides. */

  scalar fuel  = gas->YList[0];
  scalar inert = gas->YList[1];
  scalar TG = gas->T;

  fuel[top] = dirichlet (0.);
  fuel[left] = dirichlet (0.);
  fuel[right] = dirichlet (0.);

  inert[top] = dirichlet (1.);
  inert[left] = dirichlet (1.);
  inert[right] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[left] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);

  /**
  The scalar tracer is initialized with the value of the
  horizontal coordinate. */

  foreach()
    tr[] = x*f[];

  trmin = statsf(tr).min;
  trmax = statsf(tr).max;
}

/**
We adapt the grid according to the fuel mass fraction,
the temperature, and the velocity fields. */

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

event logger (t += 0.02) {
  double R = pow (3./4./pi*(2.*pi*statsf(f).sum - volumecorr), 1./3.);
  fprintf (stderr, "%d %g %.5g\n", i, t, R / (0.5*D0));
}

/**
### Grashof Number

We compute the Grashof number to quantify the intensity of the
natural convective fluxes. */

struct Grashof {
  double rhob, rhos;
  double r, g, nu;
  double value;
};

struct Grashof Gr;

event grashof (i++) {
  if (i == 0) {
    ThermoState tsg;
    tsg.T = TG0;
    tsg.P = Pref;
    tsg.x = gas->ts0->x;

    Gr.rhob = tp2.rhov (&tsg);
    Gr.nu = tp2.muv (&tsg)/tp2.rhov (&tsg);
    effective_radius = effective_radius0;
  }
  Gr.r = effective_radius;
  Gr.g = fabs (GRAVITY);

  scalar YGIntFuel = gas_int->YList[0];
  scalar YGIntInert = gas_int->YList[1];

  double TIntAvg = avg_interface (gas_int->T, f, tol=0.1);
  double YIntAvgFuel = avg_interface (YGIntFuel, f, tol=0.1);
  double YIntAvgInert = avg_interface (YGIntInert, f, tol=0.1);

  double YIntAvg[] = {YIntAvgFuel, YIntAvgInert};
  double XIntAvg[NGS];

  correctfrac (YIntAvg, NGS);
  mass2molefrac (XIntAvg, YIntAvg, gas->MWs, NGS);

  ThermoState tsg;
  tsg.P = Pref;
  if (i == 0) {
    tsg.T = TL0;
    tsg.x = gas->ts0->x;
  }
  else {
    tsg.T = TIntAvg;
    tsg.x = XIntAvg;
  }
  Gr.rhos = tp2.rhov (&tsg);

  Gr.value = (Gr.rhos - Gr.rhob)*pow (Gr.r, 3.)*Gr.g/(Gr.rhob*sq(Gr.nu));
}

/**
### Output Files

We write on a file the squared diameter decay and the
dimensionless time. */

event output_data (i += 50) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  effective_radius = pow(3./4./pi*(2.*pi*statsf(f).sum - volumecorr), 1./3.);
  double d_over_d02_old = d_over_d02;
  double tad_old = tad;

  d_over_d02 = sq (effective_radius / effective_radius0);
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
  We compute and print useful average quantities such as the
  average interface temperature and mass fractions, and the
  average droplet temperature. */

  scalar YGIntFuel = gas_int->YList[0];
  double TIntAvg = avg_interface (gas_int->T, f, tol=0.1);
  double YIntAvg = avg_interface (YGIntFuel, f, tol=0.1);

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

  fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g\n", t, tad, effective_radius0,
      effective_radius, d_over_d02, mLiq/mLiq0, kv, Gr.value, TIntAvg, YIntAvg,
      TDropAvg);
  fflush (fp);
}

/**
### Movie

We write the animation with the evolution of the
fuel mass fraction, the interface position
and the temperature field. */

#if MOVIE
event movie (t += MOVIE_EVERY) {
  clear();
  box();
  view (tx = 0.025, fov = 3.5, samples = 2,
      theta=0., phi=0., psi=-pi/2.);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = TL0, max = TG0, linear = true);
  save ("temperature.mp4");

  clear();
  box();
  view (tx = 0.025, fov = 3.5, samples = 2,
      theta=0., phi=0., psi=-pi/2.);
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
  if (i > 1) {
    char name[80];
    sprintf (name, "snapshots-%g", t);
    dump (name);
  }
}
#endif

/**
### Stopping Condition

We stop the simulation when the droplet is almost fully consumed. */

event stop (i++) {
  if (d_over_d02 <= MAX_DD02)
    return 1;
}

#if TEST_RESTART
event dump (t = 0.05) {
  dump ("restart");
  return 1;
}
#endif

/**
We run the simulation for long time, which is not reached because
the stopping condition on the droplet diameter is reached first. */

#if BASILISK_SERVER
event end (t = 0.5);
#else
event end (t = 50.);
#endif

