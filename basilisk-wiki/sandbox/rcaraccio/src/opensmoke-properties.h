/**
## OpenSMOKE++ Properties Interface
We link the OpenSMOKE++ interface to read thermodynamic and transport properties for the gas and solid phases. 
The properties are computed as a function of the local thermodynamic state, which is defined by the temperature, pressure and composition.
*/

#include "opensmoke.h"
#include "variable-properties.h"

double opensmoke_gasprop_density (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  double MWmix = OpenSMOKE_MolecularWeight_From_MoleFractions (ts->x);
  return OpenSMOKE_GasProp_Density_IdealGas (ts->T, ts->P, MWmix);
}

double opensmoke_gasprop_viscosity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_DynamicViscosity (ts->x);
}

double opensmoke_gasprop_thermalconductivity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_ThermalConductivity (ts->x);
}

double opensmoke_gasprop_heatcapacity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_HeatCapacity (ts->x);
}

void opensmoke_gasprop_diff (void * p, double * Dmix) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  OpenSMOKE_GasProp_Dmix (ts->x, Dmix);
}

void opensmoke_gasprop_species_heatcapacity (void * p, double * cpi) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  OpenSMOKE_GasProp_HeatCapacity_PureSpecies (cpi);
}

double opensmoke_solprop_heatcapacity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_SolProp_SetTemperature (ts->T);
  OpenSMOKE_SolProp_SetPressure (ts->P);
  return OpenSMOKE_SolProp_HeatCapacity (ts->x);
}

/**
### *opensmoke_solprop_density*: dummy function, rhos is constant and is not read from OS
*/

extern double rhoS;

trace
double opensmoke_solprop_density (void * p) {
  return rhoS;
}

/**
### *opensmoke_solprop_thermalconductivity()*: dummy function, lambdas is constant and is not read from OS
*/

extern double lambdaS;

double opensmoke_solprop_thermalconductivity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  
  // Corbetta wood conductivity model
  // double char_cond = 0.1405;
  // double bio_cond  = 0.1937;

  // Anca Couce wood conductivity model
  double char_cond = 0.125;
  double bio_cond = 0.056 + 2.6*1e-4*ts->T;

  double char_fraction = ts->x[OpenSMOKE_IndexOfSolidSpecies("CHAR")];
  return char_cond*char_fraction + bio_cond*(1.-char_fraction);
  // return lambdaS; //const case
}

/**
## Thermodynamic Properties

We create the instance of two structures with the
thermodynamic properties, *tp1* for the solid phase
and *tp2* for the gas phase. The same nomenclature is used
for the thermodynamic states.
*/

ThermoProps tpS, tpG;

event defaults (i = 0) {

  tpS.rhov    = opensmoke_solprop_density;
  tpS.lambdav = opensmoke_solprop_thermalconductivity;
  tpS.cpv     = opensmoke_solprop_heatcapacity;

  tpG.rhov    = opensmoke_gasprop_density;
  tpG.muv     = opensmoke_gasprop_viscosity;
  tpG.lambdav = opensmoke_gasprop_thermalconductivity;
  tpG.cpv     = opensmoke_gasprop_heatcapacity;
  tpG.diff    = opensmoke_gasprop_diff;
  tpG.cpvs    = opensmoke_gasprop_species_heatcapacity;
}
