#include "cantera/cantera.h"
#include "common-chemistry.h"

typedef void(*odefunction)(const double * y, const double t, double * dy, void * args);
void batch_isothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {}
void batch_nonisothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {}

void cantera_ode_solver (ode_function batch,
    unsigned int NEQ, double dt, double * y0, void * args)
{
  size_t ns = thermo_nSpecies (thermo);
  bool isothermal = (ns == NEQ) ? true : false;

  UserDataODE * data = (UserDataODE *)args;
  double P = data->P;
  double T = isothermal ? data->T : y0[ns];
  double * sources = data->sources;

  double ymass[ns];
  for (int i = 0; i < ns; i++)
    ymass[i] = y0[i];

  thermo_setTemperature (thermo, T);
  thermo_setPressure (thermo, P);
  thermo_setMassFractions (thermo, ns, ymass, 1);

  int reactor = reactor_new ("IdealGasConstPressureReactor", soln, "batch");
  reactor_setEnergy (reactor, !isothermal);
  int net = reactornet_new();
  reactornet_addreactor (net, reactor);
  reactornet_setTolerances (net, 1e-5, 1e-7);

  reactornet_advance (net, dt);

  // Recover results
  for (int i = 0; i < ns; i++) {
    ymass[i] = reactor_massFraction (reactor, i);
    y0[i] = ymass[i];
  }
  if (!isothermal) {
    T = reactor_temperature (reactor);
    y0[ns] = T;
  }

  // Computes derivatives (RHS) for the expansion term
  thermo_setTemperature (thermo, T);
  thermo_setPressure (thermo, P);
  thermo_setMassFractions (thermo, ns, ymass, 1);

  double wdot[ns], hm[ns], MW[ns];
  thermo_getPartialMolarEnthalpies (thermo, ns, hm);  // [J/kmol]
  kin_getNetProductionRates (kin, ns, wdot);          // [kmol/m3/s]
  thermo_getMolecularWeights (thermo, ns, MW);

  sources[ns] = 0.;
  for (int i = 0; i < ns; i++) {
    sources[i] = wdot[i]*MW[i];
    sources[ns] -= wdot[i]*hm[i];
  }

  // Re-assign values to the UserDataODE, it may be useful for binning
  data->rho = thermo_density (thermo);
  data->cp = thermo_cp_mass (thermo);

  reactor_del (reactor);
  reactornet_del (net);
}

event defaults (i = 0) {
  stiff_ode_solver = cantera_ode_solver;
}
