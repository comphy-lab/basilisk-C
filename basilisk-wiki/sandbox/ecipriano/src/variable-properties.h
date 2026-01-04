/**
# Variable Properties

Simulations involving two-phase flows with variable material
properties can be performed using this module, which defines
structures and functions that help to setup such cases.
*/

#ifndef THERMODYNAMICS_H
# define THERMODYNAMICS_H

#define VARPROP 1

#include "thermodynamics.h"

/**
## Fields Allocations

We allocate fields required by this module. All properties
are intialized as constant fields, and they are initialized
just by the solver that needs them. If a thermal solver is
used, there is no need to initialize a non-constant
diffusivity of the chemical species, for example.
*/

//(const) scalar rho1v = zeroc, rho2v = zeroc;
//(const) scalar mu1v = zeroc, mu2v = zeroc;
//(const) scalar lambda1v = zeroc, lambda2v = zeroc;
//(const) scalar cp1v = zeroc, cp2v = zeroc;

/**
## Thermodynamic State

We define a structure with variables that define the thermodynamic
state of the mixture: temperature *T*, pressure *P*, composition *x*.
*/

typedef struct {
  double T, P;
  double * x;
} ThermoState;

/**
## Thermodynamic Properties

We define a struct that contains a bunch of function pointers
for each property of interest for the gas phase.
*/

typedef struct {
  // Mixture properties
  double (* rhov) (void *);
  double (* muv) (void *);
  double (* lambdav) (void *);
  double (* cpv) (void *);
  // Species properties
  void (* dhev) (void *, double *);
  void (* diff) (void *, double *);
  void (* cps) (void *, double *);
  void (* sigmas) (void *, double *);
  // Expansion functions
  double (* betaT) (const void *, void *);
  void (* betaY) (const void *, void *, double *);
} ThermoProps;

// Functions for simpler use of ThermoState

ThermoState * new_thermo_state (size_t n) {
  ThermoState * ts = malloc (sizeof (ThermoState));
  ts->x = malloc (n*sizeof (double));
  return ts;
}

void free_thermo_state (ThermoState * ts) {
  free (ts->x), ts->x = NULL;
  free (ts), ts = NULL;
}

void copy_thermo_state (ThermoState * dest, const ThermoState * orig,
    size_t n)
{
  dest->T = orig->T;
  dest->P = orig->P;
  if (orig->x) {
    for (size_t i = 0; i < n; i++)
      dest->x[i] = orig->x[i];
  }
}

/**
## Useful functions

We define functions that are useful for variable properties
simulations.
*/

/**
### *check_termostate()*: check that the thermodynamic state is
reasonable. */

//int check_thermostate (ThermoState * ts, int NS) {
//  double sum = 0.;
//  for (int jj=0; jj<NS; jj++)
//    sum += ts->x[jj];
//
//  int T_ok = (ts->T > 180. && ts->T < 4000.) ? true : false;
//  int P_ok = (ts->P > 1e3 && ts->P < 1e7) ? true : false;
//  int X_ok = (sum > 1.-1.e-3 && sum < 1.+1.e-3) ? true : false;
//
//  return T_ok*P_ok*X_ok;
//}

/**
### *print_thermostate()*: print the thermodynamic state of the mixture
*/

void print_thermostate (ThermoState * ts, int NS, FILE * fp = stdout) {
  fprintf (fp, "Temperature = %g - Pressure = %g\n", ts->T, ts->P);
  for (int jj=0; jj<NS; jj++)
    fprintf (fp, "  Composition[%d] = %g\n", jj, ts->x[jj]);
  fprintf (fp, "\n");
}

/**
### *print_thermoprop()*: print the thermodynamic properties of the mixture
*/

void print_thermoprop (ThermoProps * tp, ThermoState * ts, int NS,
    FILE * fp = stdout)
{
  if (tp->rhov) fprintf (fp, "density = %g\n", tp->rhov (ts));
  if (tp->muv) fprintf (fp, "viscosity = %g\n", tp->muv (ts));
  if (tp->lambdav) fprintf (fp, "lambda = %g\n", tp->lambdav (ts));
  if (tp->cpv) fprintf (fp, "cp = %g\n", tp->cpv (ts));

  double dhev[NS], diff[NS], cps[NS], sigmas[NS];
  if (tp->dhev) {
    tp->dhev (ts, dhev);
    for (int i = 0; i < NS; i++)
      fprintf (fp, "dhev[%d] = %g\n", i, dhev[i]);
  }
  if (tp->diff) {
    tp->diff (ts, diff);
    for (int i = 0; i < NS; i++)
      fprintf (fp, "diff[%d] = %g\n", i, diff[i]);
  }
  if (tp->cps) {
    tp->cps (ts, cps);
    for (int i = 0; i < NS; i++)
      fprintf (fp, "cps[%d] = %g\n", i, cps[i]);
  }
  if (tp->sigmas) {
    tp->sigmas (ts, sigmas);
    for (int i = 0; i < NS; i++)
      fprintf (fp, "sigmas[%d] = %g\n", i, sigmas[i]);
  }
  fprintf (fp, "\n");
}

/**
### *gasprop_thermal_expansion()*: Thermal expansion coefficient of an ideal gas
*/

//double gasprop_thermal_expansion (ThermoProps * tp, ThermoState * ts) {
//  return ts->T > 0. ? 1./ts->T : 0.;
//}

//double gasprop_thermal_expansion (void * p, void * s) {
//  ThermoState * ts = (ThermoState *)s;
//  return ts->T > 0. ? 1./ts->T : 0.;
//}

/**
### *gasprop_species_expansion()*: Thermal expansion coefficient of an ideal gas
*/

//void gasprop_species_expansion (void * p, void * s, double * r) {
//}

/**
### *liqprop_thermal_expansion()*: Thermal expansion coefficient of a liquid
*/

//double liqprop_thermal_expansion (ThermoProps * tp, ThermoState * ts) {
//  double epsT = 1.e-3;
//  double Ttop = ts->T + epsT, Tbot = ts->T - epsT;
//  ThermoState tstop, tsbot;
//  tstop.T = Ttop, tstop.P = ts->P, tstop.x = ts->x;
//  tsbot.T = Tbot, tsbot.P = ts->P, tsbot.x = ts->x;
//  double rhotop = tp->rhov (&tstop), rhobot = tp->rhov (&tsbot);
//  double rhoval = tp->rhov (ts);
//  return (rhoval > 0.) ? -1./rhoval*(rhotop - rhobot)/(2.*epsT) : 0.;
//}

//double liqprop_thermal_expansion (void * p, void * s) {
//  ThermoProps * tp = (ThermoProps *)p;
//  ThermoState * ts = (ThermoState *)s;
//
//  double epsT = 1.e-3;
//  double Ttop = ts->T + epsT, Tbot = ts->T - epsT;
//  ThermoState tstop, tsbot;
//  tstop.T = Ttop, tstop.P = ts->P, tstop.x = ts->x;
//  tsbot.T = Tbot, tsbot.P = ts->P, tsbot.x = ts->x;
//  double rhotop = tp->rhov (&tstop), rhobot = tp->rhov (&tsbot);
//  double rhoval = tp->rhov (ts);
//  return (rhoval > 0.) ? -1./rhoval*(rhotop - rhobot)/(2.*epsT) : 0.;
//}

/**
## *mass2molefrac()*: Compute mole fractions from mass fractions

* *X*: vector filled with mole fractions
* *W*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

void mass2molefrac (double * X, const double * W, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += W[i]/MW[i];
  }
  for (int i=0; i<NS; i++) {
    X[i] = W[i]/MW[i]/(MWmix + 1.e-10);
  }
}

/**
## *mole2massfrac()*: Compute mass fractions from mole fractions

* *W*: vector filled with mole fractions
* *X*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

void mole2massfrac (double * W, const double * X, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += X[i]*MW[i];
  }
  for (int i=0; i<NS; i++) {
    W[i] = X[i]*MW[i]/(MWmix + 1.e-10);
  }
}

/**
## *mass2mw()*: Compute mixture molecular weight from mass fractions

* *W*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

double mass2mw (const double * W, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += W[i]/MW[i];
  }
  return 1./(MWmix + 1.e-10);
}

/**
## *mole2mw()*: Compute mixture molecular weight from mole fractions

* *X*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

double mole2mw (const double * X, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += X[i]*MW[i];
  }
  return MWmix;
}

/**
## *correctfrac()*: Close to 1 a vector of mass or mole fractions

* *X*: vector with mass or mole fractions
* *NS* total number of species (vector length)
*/

void correctfrac (double * X, const int NS)
{
  double sum = 0.;
  for (int i=0; i<NS; i++)
    sum += (X[i] >= 0.) ? X[i] : 0.;
  for (int i=0; i<NS; i++)
    X[i] = (X[i] >= 0.) ? X[i]/(sum + 1.e-10) : 0.;
}


#endif
