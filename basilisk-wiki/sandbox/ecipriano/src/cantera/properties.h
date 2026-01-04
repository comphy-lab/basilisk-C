/**
# Cantera Properties

We compute the material properties of a mixture using the Cantera library.
*/

#include "cantera/cantera.h"
#include "variable-properties.h"

#define VARIABLE_PROPERTIES 1

/**
## Properties Functions

Functions for the update of the density, given the thermodynamic
state.
*/

/**
### *cantera_gasprop_density()*: gas phase density according to the ideal gas low
*/

// fixme: check discrepancy between ideal gas and `thermo_density()`
double cantera_gasprop_density (void * p) {
  ThermoState * ts = (ThermoState *)p;
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  size_t ns = thermo_nSpecies (thermo);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);
  //return thermo_density (thermo);
  double MWmix = thermo_meanMolecularWeight (thermo);
  return ts->P*MWmix / (R_GAS*1e3*ts->T);
}

/**
### *cantera_gasprop_viscosity()*: gas phase dynamic viscosity
*/

double cantera_gasprop_viscosity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  size_t ns = thermo_nSpecies (thermo);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);
  return trans_viscosity (tran);
}

/**
### *cantera_gasprop_thermalconductivity()*: gas phase thermal conductivity
*/

double cantera_gasprop_thermalconductivity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  size_t ns = thermo_nSpecies (thermo);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);
  return trans_thermalConductivity (thermo);
}

/**
### *cantera_gasprop_heatcapacity()*: gas phase specific heat capacity
*/

double cantera_gasprop_heatcapacity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  size_t ns = thermo_nSpecies (thermo);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);
  return thermo_cp_mass (thermo);
}

/**
### *cantera_gasprop_heatcapacity_species()*: gas phase species heat capacity
*/

void cantera_gasprop_heatcapacity_species (void * p, double * r) {
  ThermoState * ts = (ThermoState *)p;
  size_t ns = thermo_nSpecies (thermo);
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);

  thermo_getPartialMolarCp (thermo, ns, r);

  double MW[ns];
  thermo_getMolecularWeights (thermo, ns, MW);

  for (int i = 0; i < ns; i++)
    r[i] /= MW[i];
}

/**
### *cantera_gasprop_diff()*: diffusion coefficient of a species in gas phase
*/

void cantera_gasprop_diff (void * p, double * r) {
  ThermoState * ts = (ThermoState *)p;
  size_t ns = thermo_nSpecies (thermo);
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  double xm[ns];
  for (int i = 0; i < ns; i++)
    xm[i] = max (ts->x[i], 1e-10);
  thermo_setMoleFractions (thermo, ns, xm, 1);
  trans_getMixDiffCoeffs (tran, ns, r);
}

/**
### *cantera_antoine()*: implementation of the antoine function using opensmoke
*/

double cantera_antoine (double T, double P, int i) {
  return 0.;
}

/**
### *cantera_gasprop_thermal_expansion()*: gas thermal expansion coefficient
*/

double cantera_gasprop_thermal_expansion (const void * p, void * s) {
  ThermoState * ts = (ThermoState *)s;
  return (ts->T > 0.) ? 1./ts->T : 0.;
}

/**
### *cantera_gasprop_species_expansion()*: gas species expansion coefficient
*/

void cantera_gasprop_species_expansion (const void * p, void * s, double * r) {
  ThermoState * ts = (ThermoState *)s;
  thermo_setTemperature (thermo, ts->T);
  thermo_setPressure (thermo, ts->P);
  size_t ns = thermo_nSpecies (thermo);
  thermo_setMoleFractions (thermo, ns, ts->x, 1);
  double MWmix = thermo_meanMolecularWeight (thermo);

  double MW[ns];
  thermo_getMolecularWeights (thermo, ns, MW);

  for (int i = 0; i < ns; i++)
    r[i] = MWmix / MW[i];
}

/**
### *cantera_liqprop_thermal_expansion()*: liq thermal expansion coefficient
*/

double cantera_liqprop_thermal_expansion (const void * p, void * s) {
  ThermoProps * tp = (ThermoProps *)p;
  ThermoState * ts = (ThermoState *)s;

  if (tp->rhov == NULL)
    return 0;
  else {
    double epsT = 1.e-3;
    double Ttop = ts->T + epsT, Tbot = ts->T - epsT;
    ThermoState tstop, tsbot;
    tstop.T = Ttop, tstop.P = ts->P, tstop.x = ts->x;
    tsbot.T = Tbot, tsbot.P = ts->P, tsbot.x = ts->x;
    double rhotop = tp->rhov (&tstop), rhobot = tp->rhov (&tsbot);
    double rhoval = tp->rhov (ts);
    return (rhoval > 0.) ? -1./rhoval*(rhotop - rhobot)/(2.*epsT) : 0.;
  }
}

/**
### *cantera_liqprop_species_expansion()*: liq species expansion coefficient
*/

void cantera_liqprop_species_expansion (const void * p, void * s, double * r)
{
  for (unsigned int i = 0; i < thermo_nSpecies (thermo_liq); i++)
    r[i] = 0.;
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

  tp1.rhov    = NULL;
  tp1.muv     = NULL;
  tp1.lambdav = NULL;
  tp1.cpv     = NULL;
  tp1.dhev    = NULL;
  tp1.diff    = NULL;
  tp1.cps     = NULL;
  tp1.sigmas  = NULL;
  tp1.betaT   = cantera_liqprop_thermal_expansion;
  tp1.betaY   = cantera_liqprop_species_expansion;

  tp2.rhov    = cantera_gasprop_density;
  tp2.muv     = cantera_gasprop_viscosity;
  tp2.lambdav = cantera_gasprop_thermalconductivity;
  tp2.cpv     = cantera_gasprop_heatcapacity;
  tp2.diff    = cantera_gasprop_diff;
  tp2.cps     = cantera_gasprop_heatcapacity_species;
  tp2.betaT   = cantera_gasprop_thermal_expansion;
  tp2.betaY   = cantera_gasprop_species_expansion;
}

