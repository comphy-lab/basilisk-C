/**
# Variable Properties

Simulations involving two-phase flows with variable material
properties can be performed using this module, which defines
structures and functions that help to setup such cases.
*/

#ifndef VARIABLE_PROPERTIES_H
# define VARIABLE_PROPERTIES_H

#define VARPROP 1
#include "thermodynamics.h"

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
  double (* chiT) (const void *, void *);
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

For an ideal gas the density is inversely proportional to the temperature,
therefore the thermal expansion coefficient reduces to the inverse of the
temperature, without any need for a numerical differentiation.
*/

double gasprop_thermal_expansion (const void * p, void * s) {
  ThermoState * ts = (ThermoState *)s;
  return (ts->T > 0.) ? 1./ts->T : 0.;
}

/**
### *liqprop_thermal_expansion()*: Thermal expansion coefficient of a liquid

The thermal expansion coefficient:
$$
  \beta_T = -\dfrac{1}{\rho}
  \left(\dfrac{\partial\rho}{\partial T}\right)_{P,\omega}
$$
is obtained from the numerical differentiation of the density with respect
to the temperature. The density is reached through the *rhov* function
pointer, therefore this implementation does not depend on the specific
thermodynamic backend, and the backends just forward to this function.
*/

double liqprop_thermal_expansion (const void * p, void * s) {
  ThermoProps * tp = (ThermoProps *)p;
  ThermoState * ts = (ThermoState *)s;

  if (tp->rhov == NULL)
    return 0.;
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
## Isothermal Compressibility

The isothermal compressibility:
$$
  \chi_T = \dfrac{1}{\rho}\left(\dfrac{\partial\rho}{\partial P}\right)_{T,\omega}
$$
is required by closed systems, where the thermodynamic pressure changes in time
and it contributes to the divergence of the velocity field. It is the pressure
counterpart of the thermal expansion coefficient *betaT*.
*/

/**
### *gasprop_isothermal_compressibility()*: Isothermal compressibility of an ideal gas
*/

double gasprop_isothermal_compressibility (const void * p, void * s) {
  ThermoState * ts = (ThermoState *)s;
  return (ts->P > 0.) ? 1./ts->P : 0.;
}

/**
### *liqprop_isothermal_compressibility()*: Isothermal compressibility of a liquid

The compressibility of a generic phase is obtained from the numerical
differentiation of the density with respect to the pressure. The perturbation
is relative, because the absolute value of the pressure can be large.
*/

double liqprop_isothermal_compressibility (const void * p, void * s) {
  ThermoProps * tp = (ThermoProps *)p;
  ThermoState * ts = (ThermoState *)s;

  if (tp->rhov == NULL)
    return 0.;
  else {
    double epsP = 1.e-4*ts->P;
    ThermoState tstop, tsbot;
    tstop.T = ts->T, tstop.P = ts->P + epsP, tstop.x = ts->x;
    tsbot.T = ts->T, tsbot.P = ts->P - epsP, tsbot.x = ts->x;
    double rhotop = tp->rhov (&tstop), rhobot = tp->rhov (&tsbot);
    double rhoval = tp->rhov (ts);
    return (rhoval > 0. && epsP > 0.) ?
      1./rhoval*(rhotop - rhobot)/(2.*epsP) : 0.;
  }
}
#endif
