/**
# Basilisk Properties

We define functions for the calculation of variable properties. These
functions are specific to the chemical species being solved, and they
need to be assigned to the `ThermoProp` pointers depending on the
problem.

This file is used to test the variable properties formulation on the
Basilisk server, and it contains functions which are tailored for the
problem being solved. For example, it implements the functions for a
specific chemical species or mixture, and it is unable to manage any
possible mixture, which would require a proper library (see for example
[opensmoke-properties.h](/sandbox/ecipriano/src/opensmoke-properties.h).
*/

#include "thermodynamics.h"
#include "fractions.h"
#include "common-evaporation.h"
#include "variable-properties.h"

/**
## Properties Functions

Functions for the update of the density, given the thermodynamic
state.
*/

extern double inMW[NGS];

/**
## Constant properties

We set functions that returns the constant properties, using the same
variables already in use for the constant properties simulations. */

double const_liqprop_density (void * p) {
  extern double rho1;
  return rho1;
}

double const_liqprop_viscosity (void * p) {
  extern double mu1;
  return mu1;
}

double const_liqprop_thermalconductivity (void * p) {
  extern double lambda1;
  return lambda1;
}

double const_liqprop_heatcapacity (void * p) {
  extern double cp1;
  return cp1;
}

double const_liqprop_heatcapacity_species (void * p, int i) {
  extern double cp1;
  return cp1;
}

double const_liqprop_dhev (void * p, int i) {
  extern double dhev;
  return dhev;
}

double const_liqprop_diff (void * p, int i) {
  extern double inDmix1[NLS];
  return inDmix1[i];
}

double const_gasprop_density (void * p) {
  extern double rho2;
  return rho2;
}

double const_gasprop_viscosity (void * p) {
  extern double mu2;
  return mu2;
}

double const_gasprop_thermalconductivity (void * p) {
  extern double lambda2;
  return lambda2;
}

double const_gasprop_heatcapacity (void * p) {
  extern double cp2;
  return cp2;
}

double const_gasprop_heatcapacity_species (void * p, int i) {
  extern double cp2;
  return cp2;
}

double const_gasprop_diff (void * p, int i) {
  extern double inDmix2[NGS];
  return inDmix2[i];
}

/**
## Variable Properties Functions

Other functions implementing laws for variable properties. */

double gasprop_density_idealgas (void * p) {
  ThermoState * ts = (ThermoState *)p;
  double MWmix = mole2mw (ts->x, inMW, NGS);
  return ts->P*MWmix/(R_GAS*1e3)/ts->T;
}

double gasprop_viscosity_heptane (void * p) {
  double c[] = {-2.529079761891171074e+01, 3.936332538558730665e+00, -3.545774416702855980e-01, 1.250482423406749401e-02};
  ThermoState * ts = (ThermoState *)p;
  const double logT = log (ts->T);
  const double logT2 = logT*logT;
  const double logT3 = logT*logT2;
  return exp (c[0] + c[1]*logT + c[2]*logT2 + c[3]*logT3);
}

double gasprop_viscosity_nitrogen (void * p) {
  double c[] = {-1.853769421193775102e+01, 2.232970914651082772e+00, -2.102014061222551300e-01, 9.264621890951153793e-03};
  ThermoState * ts = (ThermoState *)p;
  const double logT = log (ts->T);
  const double logT2 = logT*logT;
  const double logT3 = logT*logT2;
  return exp (c[0] + c[1]*logT + c[2]*logT2 + c[3]*logT3);
}

double gasprop_viscosity_heptane_nitrogen (void * p) {
  double eta[2] = {gasprop_viscosity_heptane (p), gasprop_viscosity_nitrogen (p)};
  ThermoState * ts = (ThermoState *)p;

  double sum = 0;
  for (unsigned int k = 0; k < NGS; k++)
    sum += ts->x[k] * sqrt (inMW[k]);

  double etamix = 0.;
  for (unsigned int k = 0; k < NGS; k++)
    etamix += ts->x[k] * eta[k] * sqrt (inMW[k]);
  etamix /= sum;

  return etamix;
}

double liqprop_density_heptane (void * p) {
  double a = 5.2745973, b = 0.07741, c = 557.342, d = 0.13673;
  ThermoState * ts = (ThermoState *)p;
  double Teff = min (ts->T, 0.999*c);
  return a / (pow (b, 1. + pow (1. - Teff / c, d)));
}

double liqprop_viscosity_heptane (void * p) {
  double a = -12.805, b = 1237.54, c = 0.001869, d = 0.002546e-3, e = 0.;
  ThermoState * ts = (ThermoState *)p;
  return exp (a + b/ts->T + c*ts->T + d*pow (ts->T, 2.) + e*pow (ts->T, 3.));
}

double liqprop_thermalconductivity_heptane (void * p) {
  double a = 0.2149, b = -0.02998e-02, c = -0.3e-07, d = 0.103e-09, e = -0.12e-12;
  double Tc = 540.3;
  ThermoState * ts = (ThermoState *)p;
  double Teff = min (ts->T, 0.999*Tc);
  return (((((((e * Teff) + d) * Teff) + c) * Teff) + b) * Teff + a);
}

double liqprop_heatcapacity_heptane_species (void * p, int i) {
  double a = 1925.0, b = -0.238, c = 0.001105, d = 0.1081e-04, e = 0.;
  ThermoState * ts = (ThermoState *)p;
  return a + b*ts->T + c*pow (ts->T, 2.) + d*pow (ts->T, 3.) + e*pow (ts->T, -2.);
}

double liqprop_dhev_heptane (void * p, int i) {
  double a = 499720.0, b = 0.39266, c = -0.005169, d = -0.003931, e = 0.005262;
  double Tc = 540.3;
  ThermoState * ts = (ThermoState *)p;
  double TTc = min (ts->T, Tc)/Tc;
  return a*pow ((1. - TTc), b + c*TTc + d*pow (TTc, 2.) + e*pow (TTc, 3.));
}

double liqprop_sigma_heptane (void * p, int i) {
  double a = 53.64, b = 1.2431;
  double Tc = 540.3;
  ThermoState * ts = (ThermoState *)p;
  return 1.e-3*a*pow (1. - min (ts->T, Tc)/Tc, b);
}

/**
## Thermodynamic Properties

We create the instance of two structures with the
thermodynamic properties, *tp1* for the liquid phase
and *tp2* for the gas phase. The same nomenclature is used
for the thermodynamic states.
*/

ThermoProps tp1, tp2;

/**
## Initialization

We set the thermodynamic properties function pointers
to the specific opensmoke functions declared above.
*/

event defaults (i = 0) {

  /**
  We set the thermodynamic properties functions to the
  correct opensmoke functions that compute material
  properties. */

  tp1.rhov    = const_liqprop_density;
  tp1.muv     = const_liqprop_viscosity;
  tp1.lambdav = const_liqprop_thermalconductivity;
  tp1.cpv     = const_liqprop_heatcapacity;
  tp1.pvap    = NULL;
  tp1.dhev    = const_liqprop_dhev;
  tp1.diff    = const_liqprop_diff;
  tp1.cps     = const_liqprop_heatcapacity_species;
  tp1.sigmas  = NULL;

  tp2.rhov    = const_gasprop_density;
  tp2.muv     = const_gasprop_viscosity;
  tp2.lambdav = const_gasprop_thermalconductivity;
  tp2.cpv     = const_gasprop_heatcapacity;
  tp2.diff    = const_gasprop_diff;
  tp2.cps     = const_gasprop_heatcapacity_species;
}

