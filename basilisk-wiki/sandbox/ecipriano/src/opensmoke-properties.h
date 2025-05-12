/**
# OpenSMOKE++ Properties

We compute the material properties of a mixture using the
OpenSMOKE++ library.
*/

#include "opensmoke.h"
#include "variable-properties.h"

/**
## Properties Functions

Functions for the update of the density, given the thermodynamic
state.
*/

/**
### *opensmoke_gasprop_density()*: gas phase density according to the ideal gas low
*/

double opensmoke_gasprop_density (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  double MWmix = OpenSMOKE_MolecularWeight_From_MoleFractions (ts->x);
  return OpenSMOKE_GasProp_Density_IdealGas (ts->T, ts->P, MWmix);
}

/**
### *opensmoke_gasprop_viscosity()*: gas phase dynamic viscosity
*/

double opensmoke_gasprop_viscosity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_DynamicViscosity (ts->x);
}

/**
### *opensmoke_gasprop_thermalconductivity()*: gas phase thermal conductivity
*/

double opensmoke_gasprop_thermalconductivity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_ThermalConductivity (ts->x);
}

/**
### *opensmoke_gasprop_heatcapacity()*: gas phase specific heat capacity
*/

double opensmoke_gasprop_heatcapacity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_HeatCapacity (ts->x);
}

/**
### *opensmoke_gasprop_heatcapacity_species()*: gas phase species heat capacity
*/

double opensmoke_gasprop_heatcapacity_species (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  double Cps[OpenSMOKE_NumberOfSpecies()];
  OpenSMOKE_GasProp_HeatCapacity_PureSpecies (Cps);
  return Cps[i];
}

/**
### *opensmoke_gasprop_diff()*: diffusion coefficient of a species in gas phase
*/

double opensmoke_gasprop_diff (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_Dmix (ts->x, i);
}

/**
### *opensmoke_liqprop_density_addvol()*: liquid phase mixture density with the additive volume method
*/

double opensmoke_liqprop_density_addvol (void * p) {
  ThermoState * ts = (ThermoState *)p;
  return OpenSMOKE_LiqProp_Density_Mix_AddVol (ts->T, ts->P, ts->x);
}

/**
### *opensmoke_liqprop_viscosity()*: liquid phase mixture dynamic viscosity
*/

double opensmoke_liqprop_viscosity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  return OpenSMOKE_LiqProp_DynamicViscosity_Mix (ts->T, ts->x);
}

/**
### *opensmoke_liqprop_thermalconductivity()*: liquid phase mixture thermal conductivity
*/

double opensmoke_liqprop_thermalconductivity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  return OpenSMOKE_LiqProp_ThermalConductivity_Mix (ts->T, ts->x);
}

/**
### *opensmoke_liqprop_heatcapacity()*: liquid phase mixture specific heat capacity
*/

double opensmoke_liqprop_heatcapacity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  return OpenSMOKE_LiqProp_HeatCapacity_Mix (ts->T, ts->x);
}

/**
### *opensmoke_liqprop_heatcapacity_species()*: liquid phase species heat capacity
*/

double opensmoke_liqprop_heatcapacity_species (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  const char* name = OpenSMOKE_NamesOfLiquidSpecies (i);
  int len = strlen (name);
  char corrname[len+1];
  strcpy (corrname, name);
  corrname[3 <= len ? len-3 : 0] = '\0';
  return OpenSMOKE_LiqProp_HeatCapacity_PureSpecies (corrname, ts->T);
}

/**
### *opensmoke_liqprop_pvap()*: vapor pressure of the chemical species
*/

double opensmoke_liqprop_pvap (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  const char* name = OpenSMOKE_NamesOfLiquidSpecies (i);
  int len = strlen (name);
  char corrname[len+1];
  strcpy (corrname, name);
  corrname[3 <= len ? len-3 : 0] = '\0';
  return OpenSMOKE_LiqProp_VaporPressure (corrname, ts->T, ts->P);
}

/**
### *opensmoke_liqprop_dhev()*: vapor pressure of the chemical species
*/

double opensmoke_liqprop_dhev (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  const char* name = OpenSMOKE_NamesOfLiquidSpecies (i);
  int len = strlen (name);
  char corrname[len+1];
  strcpy (corrname, name);
  corrname[3 <= len ? len-3 : 0] = '\0';
  return OpenSMOKE_LiqProp_VaporizationEnthalpy (corrname, ts->T);
}

/**
### *opensmoke_liqprop_sigma()*: surface tension of the chemical species
*/

double opensmoke_liqprop_sigma (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  const char* name = OpenSMOKE_NamesOfLiquidSpecies (i);
  int len = strlen (name);
  char corrname[len+1];
  strcpy (corrname, name);
  corrname[3 <= len ? len-3 : 0] = '\0';
  return OpenSMOKE_LiqProp_Sigma (corrname, ts->T);
}

/**
### *opensmoke_liqprop_diff_pg()*: diffusion coefficient of a species in liquid phase (Perkins Geankopolis)
*/

double opensmoke_liqprop_diff_pg (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_LiqProp_Dmix_PerkinsGeankopolis (ts->T, ts->P, ts->x, i);
}

/**
### *opensmoke_liqprop_diff_c()*: diffusion coefficient of a species in liquid phase (Cullinan)
*/

double opensmoke_liqprop_diff_c (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_LiqProp_Dmix_Cullinan (ts->T, ts->P, ts->x, i);
}

/**
### *opensmoke_liqprop_diff_lc()*: diffusion coefficient of a species in liquid phase (Leffler Cullinan)
*/

double opensmoke_liqprop_diff_lc (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_LiqProp_Dmix_LefflerCullinan (ts->T, ts->P, ts->x, i);
}

/**
### *opensmoke_antoine()*: implementation of the antoine function using opensmoke
*/

double opensmoke_antoine (double T, double P, int i) {
  ThermoState ts;
  ts.T = T, ts.P = P;
  return opensmoke_liqprop_pvap (&ts, i)/P;
}

/**
## Thermodynamic Properties

We create the instance of two structures with the
thermodynamic properties, *tp1* for the liquid phase
and *tp2* for the gas phase. The same nomenclature is used
for the thermodynamic states.
*/

ThermoProps tp1, tp2;
ThermoState ts1, ts2;

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

  tp1.rhov    = opensmoke_liqprop_density_addvol;
  tp1.muv     = opensmoke_liqprop_viscosity;
  tp1.lambdav = opensmoke_liqprop_thermalconductivity;
  tp1.cpv     = opensmoke_liqprop_heatcapacity;
  tp1.pvap    = opensmoke_liqprop_pvap;
  tp1.dhev    = opensmoke_liqprop_dhev;
  tp1.diff    = opensmoke_liqprop_diff_lc;
  tp1.cps     = opensmoke_liqprop_heatcapacity_species;
  tp1.sigmas  = opensmoke_liqprop_sigma;

  tp2.rhov    = opensmoke_gasprop_density;
  tp2.muv     = opensmoke_gasprop_viscosity;
  tp2.lambdav = opensmoke_gasprop_thermalconductivity;
  tp2.cpv     = opensmoke_gasprop_heatcapacity;
  tp2.diff    = opensmoke_gasprop_diff;
  tp2.cps     = opensmoke_gasprop_heatcapacity_species;
}

