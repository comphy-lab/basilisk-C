/**
# OpenSMOKE++ Properties

We compute the material properties of a mixture using the
OpenSMOKE++ library.
*/

#include "variable-properties.h"
#define BASILISK_PROPERTIES 1

extern int NLS, NGS;

/**
## Properties Functions

Functions for the update of the density, given the thermodynamic
state.
*/

/**
### *const_gasprop_density()*: gas phase density according to the ideal gas low
*/

double const_gasprop_density (void * p) {
  extern double rho2;
  return rho2;
}

/**
### *const_gasprop_viscosity()*: gas phase dynamic viscosity
*/

double const_gasprop_viscosity (void * p) {
  extern double mu2;
  return mu2;
}

/**
### *const_gasprop_thermalconductivity()*: gas phase thermal conductivity
*/

double const_gasprop_thermalconductivity (void * p) {
  extern double lambda2;
  return lambda2;
}

/**
### *const_gasprop_heatcapacity()*: gas phase specific heat capacity
*/

double const_gasprop_heatcapacity (void * p) {
  extern double cp2;
  return cp2;
}

/**
### *const_gasprop_heatcapacity_species()*: gas phase species heat capacity
*/

void const_gasprop_heatcapacity_species (void * p, double * r) {
  extern double cp2;
  for (int i = 0; i < NGS; i++)
    r[i] = cp2;
}

/**
### *const_gasprop_diff()*: diffusion coefficient of a species in gas phase
*/

void const_gasprop_diff (void * p, double * r) {
  extern double Dmix2;
  for (int i = 0; i < NGS; i++)
    r[i] = Dmix2;
}

/**
### *const_liqprop_density_addvol()*: liquid phase mixture density with the additive volume method
*/

double const_liqprop_density_addvol (void * p) {
  extern double rho1;
  return rho1;
}

/**
### *const_liqprop_viscosity()*: liquid phase mixture dynamic viscosity
*/

double const_liqprop_viscosity (void * p) {
  extern double mu1;
  return mu1;
}

/**
### *const_liqprop_thermalconductivity()*: liquid phase mixture thermal conductivity
*/

double const_liqprop_thermalconductivity (void * p) {
  extern double lambda1;
  return lambda1;
}

/**
### *const_liqprop_heatcapacity()*: liquid phase mixture specific heat capacity
*/

double const_liqprop_heatcapacity (void * p) {
  extern double cp1;
  return cp1;
}

/**
### *const_liqprop_heatcapacity_species()*: liquid phase species heat capacity
*/

void const_liqprop_heatcapacity_species (void * p, double * r) {
  extern double cp1;
  for (int i = 0; i < NLS; i++)
    r[i] = cp1;
}

/**
### *const_liqprop_dhev()*: vapor pressure of the chemical species
*/

void const_liqprop_dhev (void * p, double * r) {
  extern double dhev;
  for (int i = 0; i < NLS; i++)
    r[i] = dhev;
}

/**
### *const_liqprop_sigma()*: surface tension of the chemical species
*/

void const_liqprop_sigma (void * p, double * r) {
  return;
}

/**
### *const_liqprop_diff()*: diffusion coefficient of a species in liquid phase
*/

void const_liqprop_diff (void * p, double * r) {
  extern double Dmix1;
  for (int i = 0; i < NLS; i++)
    r[i] = Dmix1;
}

/**
### *const_liqprop_pvap()*: vapor pressure of the chemical species
*/

double const_liqprop_pvap (void * p, int i) {
  return 0;
}

/**
### *const_antoine()*: implementation of the antoine function using opensmoke
*/

double const_antoine (double T, double P, int i) {
  return 0;
}

/**
### *const_gasprop_thermal_expansion()*: gas thermal expansion coefficient
*/

double const_gasprop_thermal_expansion (const void * p, void * s) {
  return 0;
}

/**
### *const_gasprop_species_expansion()*: gas species expansion coefficient
*/

void const_gasprop_species_expansion (const void * p, void * s, double * r) {
  for (unsigned int i = 0; i < NGS; i++)
    r[i] = 0;
}

/**
### *const_liqprop_thermal_expansion()*: liq thermal expansion coefficient
*/

double const_liqprop_thermal_expansion (const void * p, void * s) {
  ThermoProps * tp = (ThermoProps *)p;
  ThermoState * ts = (ThermoState *)s;

  double epsT = 1.e-3;
  double Ttop = ts->T + epsT, Tbot = ts->T - epsT;
  ThermoState tstop, tsbot;
  tstop.T = Ttop, tstop.P = ts->P, tstop.x = ts->x;
  tsbot.T = Tbot, tsbot.P = ts->P, tsbot.x = ts->x;
  double rhotop = tp->rhov (&tstop), rhobot = tp->rhov (&tsbot);
  double rhoval = tp->rhov (ts);
  return (rhoval > 0.) ? -1./rhoval*(rhotop - rhobot)/(2.*epsT) : 0.;
}

/**
### *const_liqprop_species_expansion()*: liq species expansion coefficient
*/

void const_liqprop_species_expansion (const void * p, void * s, double * r)
{
  for (unsigned int i = 0; i < NLS; i++)
    r[i] = 0.;  // fixme: empty
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

  tp1.rhov    = const_liqprop_density_addvol;
  tp1.muv     = const_liqprop_viscosity;
  tp1.lambdav = const_liqprop_thermalconductivity;
  tp1.cpv     = const_liqprop_heatcapacity;
  tp1.dhev    = const_liqprop_dhev;
  tp1.diff    = const_liqprop_diff;
  tp1.cps     = const_liqprop_heatcapacity_species;
  tp1.sigmas  = const_liqprop_sigma;
  tp1.betaT   = const_liqprop_thermal_expansion;
  tp1.betaY   = const_liqprop_species_expansion;

  tp2.rhov    = const_gasprop_density;
  tp2.muv     = const_gasprop_viscosity;
  tp2.lambdav = const_gasprop_thermalconductivity;
  tp2.cpv     = const_gasprop_heatcapacity;
  tp2.diff    = const_gasprop_diff;
  tp2.cps     = const_gasprop_heatcapacity_species;
  tp2.betaT   = const_gasprop_thermal_expansion;
  tp2.betaY   = const_gasprop_species_expansion;
}
