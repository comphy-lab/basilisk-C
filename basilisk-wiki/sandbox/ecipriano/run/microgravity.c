/**
# Evaporation of a Pure Droplet in Microgravity Conditions

This module reports the simulation setup for a generic pure droplet
in reduced (zero) gravity conditions. This configuration has been
used in several experiental works, in order to isolate the droplet
evaporation process from buoyancy-driven flows. Ideally, the evaporation
of the droplet is perfectly spherically-symmetric, and data from
experimental investigation can be compared with simplified 1D models.

This simulation includes variable thermodynamic and transport properties,
and the interface radiation contribution. We neglect the heat transfer
from the suspending solid fiber, which may affect experimental data if
the fiber is large. The simulation is performed by simulating a droplet
at the corner of a square domain, exploting the axial-symmetry.
*/

/**
## Default Simulation Data

The simulation setup is independent from the liquid under investigation,
and we may want to perform the same simulation studying the sensitivity
of the numerical results to the operative conditions (temperature,
pressure, ecc...) and to the liquid fuel or the ambient conditions.
To avoid code duplication, we set a bunch of default compiler variables,
which are overwritten in the `Makefile` to create different cases with
different operative conditions.

The default properties are:

* initial ambient temperature: `TEMPERATURE = 773 K`
* initial droplet temperature: `TEMPERATURE_DROPLET = 300 K`
* constant thermodynamic pressure: `PRESSURE = 10 atm`
* initial droplet diameter: `DIAMETER = 1 mm`
* emissivity of the liquid fuel: `RADIATION_INTERFACE = 0.93`
* name of the liquid fuel: `FUEL = n-heptane`
* name of the inert/ambient species: `INERT = nitrogen`
* path of the kinetics folder: `KINFOLDER = evaporation/n-heptane-in-nitrogen`
* path of the liquid properties folder: `LIQFOLDER = LiquidProperties`

these properties describe the evaporation of a n-heptane droplet in
nitrogen in microgravity conditions, including the interface radiation.
These properties can be easily changed by overwriting the compilation
variables, for example by setting `-DTEMPERATURE=473 -DPRESSURE=1` to
change the ambient temperature and pressure. The interface radiation is
suppressed by setting a null value of emissivity: `-DRADIATION_INTERFACE=0`.

Eventually, one may consider the possibility of compiling a code which
reads the properties directly from an input file (using libconfig for example).
*/

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifndef TEMPERATURE
# define TEMPERATURE 773.
#endif

#ifndef TEMPERATURE_DROPLET
# define TEMPERATURE_DROPLET 300.
#endif

#ifndef PRESSURE
# define PRESSURE 10.
#endif

#ifndef DIAMETER
# define DIAMETER 1.e-3
#endif

#ifndef FUEL
# define FUEL NC7H16
#endif

#ifndef INERT
# define INERT N2
#endif

#ifndef KINFOLDER
# define KINFOLDER evaporation/n-heptane-in-nitrogen
#endif

#ifndef LIQFOLDER
# define LIQFOLDER LiquidProperties
#endif

#ifndef RADIATION_INTERFACE
# define RADIATION_INTERFACE 0.93
#endif

/**
## Phase Change Setup

We define the names of the gas and the liquid species, and their
intial composition (in mass fractions). */

#define NGS 2
#define NLS 1

char* gas_species[NGS] = {TOSTRING(FUEL), TOSTRING(INERT)};
char* liq_species[NLS] = {TOSTRING(FUEL)};
char* inert_species[1] = {TOSTRING(INERT)};
double gas_start[NGS] = {0., 1.};
double liq_start[NLS] = {1.};

/**
We declare all the properties necessary for the multicomponent phase change model. The properties are set to null values because they are overwritten by the variable-properties formulation, which computes all the physical properties as a function of the thermodynamic state of the mixture. */

double inDmix1[NLS] = {0.};
double inDmix2[NGS] = {0.};
double inKeq[NLS] = {0.};
double lambda1 = 0.;
double lambda2 = 0.;
double dhev = 0.;
double cp1 = 0.;
double cp2 = 0.;

/**
We set the initial temperature of the liquid and of the gas phase. */

double TL0 = TEMPERATURE_DROPLET;
double TG0 = TEMPERATURE;

/**
We solve both mass fractions and temperature fields, also in the interface
jump condition. The GSL library is used for root finding operations, whose
tolerance is controlled by the variable `FSOLVE_ABSTOL`. The thermodynamic
equlibrium is computed from Antoine's law, and we solve the momentum equation
in non-conservative form by setting the variable `NO_ADVECTION_DIV` to 1.

Additional compiler variables are used to activate terms in the governing
equations. In particular: `FICK_CORRECTED` is used to force the diffusive
fluxes to close to zero, even if the diffusivity of every chemical species
is different; `MOLAR_DIFFUSION` is used to correct Fick's law considering
that the diffusivity values are mole fractions-based; `MASS_DIFFUSION_ENTHALPY`
includes the species diffusion contribution in the temperature equation. All
of them should be used.
*/

#define SOLVE_TEMPERATURE
#define USE_GSL 0
#define FSOLVE_ABSTOL 1.e-3
#define USE_ANTOINE_OPENSMOKE
#define FICK_CORRECTED
#define MOLAR_DIFFUSION
#define MASS_DIFFUSION_ENTHALPY
#define NO_ADVECTION_DIV 1

/**
## Simulation Setup

In this simulation, the Navier-Stokes equations can be solved both using the
centered solver with divergence source term, or using the centered solver with
velocity jump. The latter should be preferred since it is able to limit oscillations in the velocity field. The [two-phase.h](/src/two-phase.h)
solver is extended to variable physical properties by including a policy
for the calculation of such properties. In this case, we use correlations implemented in OpenSMOKE++ libraries.
The solution of the jump conditions and of the temperature and mass fraction
fields is performed by the multicomponent phase change model.
*/
  
#include "axi.h"
#if JUMP
# include "navier-stokes/velocity-jump.h"
#else
# include "navier-stokes/centered-evaporation.h"
# include "navier-stokes/centered-doubled.h"
#endif
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "tension.h"
#include "recoil.h"
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
the initial droplet diameter, and additional data for
post-processing. */

int maxlevel, minlevel = 2;
double D0 = DIAMETER, effective_radius0, d_over_d02 = 1., tad = 0.;
double volume0 = 0.;

int main (void) {
  
  /**
  We set the kinetics folder, which defines the species
  of the simulation, and it is used by OpenSMOKE++ for the
  calculation of the thermodynamic and transport properties.
  The path is relative to `OpenSMOKEppInterface/kinetics`.
  We do the same for the folder gathering the properties of
  the liquid species. */

  kinfolder = TOSTRING(KINFOLDER);
  liqfolder = TOSTRING(LIQFOLDER);

  /**
  We set additional simulation properties. Be careful, these are "dummy"
  properties which are overwritten by the variable-properties formulation
  already at the first iteration. To start the simulation without any
  problem of division by zero, the viscosity value should be non-null. */

  rho1 = 1.; rho2 = 1.;
  mu1 = 1.; mu2 = 1.;
  
  /**
  We set the thermodynamic pressure in SI units (Pa). */
  
  Pref = PRESSURE*101325.;

  /**
  We change the dimensions of the domain as a function of the initial
  diameter of the droplet. The domain is large with respect to the
  droplet diameter, in order to avoid the influence of the domain
  boundaries on the droplet evaporation dynamics. */

  double RR = 7.986462e+01;
  L0 = 0.5*RR*D0;

  /**
  We use a constant surface tension value. */

  f.sigma = 0.03;

  /**
  We run the simulation at different maximum levels of refinement. */

  for (maxlevel = 10; maxlevel <= 10; maxlevel++) {
    init_grid (1 << 9);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field according to `D0`. */

event init (i = 0) {
  refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);
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

  for (int jj=0; jj<NGS; jj++)
    inMW[jj] = OpenSMOKE_MW (jj);

  /**
  On the top and right sides we set Dirichlet boundary conditions
  for the temperature and mass fraction fields. */

  scalar fuel  = YGList[OpenSMOKE_IndexOfSpecies (TOSTRING(FUEL))];
  scalar inert = YGList[OpenSMOKE_IndexOfSpecies (TOSTRING(INERT))];

  fuel[top] = dirichlet (0.);
  fuel[right] = dirichlet (0.);

  inert[top] = dirichlet (1.);
  inert[right] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

/**
We adapt the grid according to the fuel mass fraction,
the temperature, and the velocity fields. */

#if TREE
event adapt (i++) {
  scalar fuel = YList[OpenSMOKE_IndexOfSpecies (TOSTRING(FUEL))];
  adapt_wavelet_leave_interface ({fuel,T,u.x,u.y}, {f},
      (double[]){1.e-2,1.e-1,1.e-1,1.e-1}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

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
  tad = t/sq(D0*1e3);

  /**
  The vaporization rate is computed according to the formula
  in Liu & Avedisian, 2011, pag. 777 bottom. */

  double kv = 0.;
  if (i > 1)
    kv = fabs ((d_over_d02 - d_over_d02_old)/(tad - tad_old));

  double mLiq = 0.;
  foreach(reduction(+:mLiq))
    mLiq += rho1v[]*f[]*dv();

  /**
  We compute and print useful average quantities such as the
  average interface temperature and mass fractions, and the
  average droplet temperature. */

  scalar YGIntFuel = YGIntList[0];
  double TIntAvg = avg_interface (TInt, f, tol=0.1);
  double YIntAvg = avg_interface (YGIntFuel, f, tol=0.1);

  int counter = 0;
  double TDropAvg = 0.;
  foreach(reduction(+:TDropAvg) reduction(+:counter)) {
    if (f[] > 1.-F_ERR) {
      counter++;
      TDropAvg += TL[];
    }
  }
  TDropAvg = (counter > 0.) ? TDropAvg/counter : 0.;

  fprintf (fp, "%g %g %g %g %g %g %g %g %g\n", t, tad, effective_radius,
      d_over_d02, mLiq/mLiq0, kv, TIntAvg, YIntAvg, TDropAvg);
}

/**
### Movie

We write the animation with the evolution of the
fuel mass fraction, the interface position
and the temperature field. */

#if MOVIE
event movie (t += 0.01) {
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
  squares (TOSTRING(FUEL), min = 0., max = 1., linear = true);
  save ("fuel.mp4");
}
#endif

/**
### Snapshots

Output dump files for restore or post-processing. */

#if DUMP
event snapshots (t += 0.1) {
  char name[80];
  sprintf (name, "snapshots-%g", t);
  dump (name);
}
#endif

/**
### Stopping Condition

We stop the simulation when the droplet is almost fully consumed. */

event stop (i++) {
  if (d_over_d02 <= 0.05)
    return 1;
}

/**
We run the simulation for long time, which is not reached because
the stopping condition on the droplet diameter is reached first. */

event end (t = 50.);
