/**
# Variable Properties

Simulations involving two-phase flows with variable material
properties can be performed using this module, which defines
structures and functions that help to setup such cases.
*/

#define VARPROP

/**
## Fields Allocations

We allocate fields required by this module. All properties
are intialized as constant fields, and they are initialized
just by the solver that needs them. If a thermal solver is
used, there is no need to initialize a non-constant
diffusivity of the chemical species, for example.
*/

(const) scalar rho1v = zeroc, rho2v = zeroc;
(const) scalar mu1v = zeroc, mu2v = zeroc;
(const) scalar lambda1v = zeroc, lambda2v = zeroc;
(const) scalar cp1v = zeroc, cp2v = zeroc;

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
  double (* pvap) (void *, int);
  double (* dhev) (void *, int);
  double (* diff) (void *, int);
  double (* cps)  (void *, int);
  double (* sigmas) (void *, int);
} ThermoProps;

/**
## Overwrite Two-Phase Properties

We define macros used to overwrite the calculation of the variable density and
viscosity fields in two-phase simulations using any method (VOF, LS, CLSVOF).
*/

#define aavg(f,v1,v2) (clamp(f,0.,1.)*(v1 - v2) + v2)
#define havg(f,v1,v2) (1./(clamp(f,0,1)*(1./(v1) - 1./(v2)) + 1./(v2)))

/**
## Update One-Field Properties

We overwrite the event `properties` to re-calculate the density and viscosity
fields after the calculation in [two-phase.h](/src/two-phase-generic.h). */

extern scalar f;
extern face vector alphav;
extern scalar rhov;
#ifdef FILTERED
extern scalar sf;
#else
# define sf f
#endif

event properties (i++) {
  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;

    double rho1vh = 0.5*(rho1v[] + rho1v[-1]);
    double rho2vh = 0.5*(rho2v[] + rho2v[-1]);

    alphav.x[] = fm.x[]/aavg (ff, rho1vh, rho2vh);

    {
      face vector muv = mu;
      double mu1vh = 0.5*(mu1v[] + mu1v[-1]);
      double mu2vh = 0.5*(mu2v[] + mu2v[-1]);

      muv.x[] = fm.x[]*aavg (ff, mu1vh, mu2vh);
    }
  }

  foreach() {
    double rho1vh = rho1v[];
    double rho2vh = rho2v[];

    rhov[] = cm[]*aavg (sf[], rho1vh, rho2vh);
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

int check_thermostate (ThermoState * ts, int NS) {
  double sum = 0.;
  for (int jj=0; jj<NS; jj++)
    sum += ts->x[jj];

  int T_ok = (ts->T > 180. && ts->T < 4000.) ? true : false;
  int P_ok = (ts->P > 1e3 && ts->P < 1e7) ? true : false;
  int X_ok = (sum > 1.-1.e-3 && sum < 1.+1.e-3) ? true : false;

  return T_ok*P_ok*X_ok;
}

/**
### *print_thermostate()*: print the thermodynamic state of the mixture.
*/

void print_thermostate (ThermoState * ts, int NS, FILE * fp = stdout) {
  fprintf (fp, "Temperature = %g - Pressure = %g\n", ts->T, ts->P);
  for (int jj=0; jj<NS; jj++)
    fprintf (fp, "  Composition[%d] = %g\n", jj, ts->x[jj]);
  fprintf (fp, "\n");
}

/**
### *print_thermoprop()*: print the thermodynamic properties of the mixture.
*/

void print_thermoprop (ThermoProps * tp, ThermoState * ts, int NS, FILE * fp = stdout) {
  if (tp->rhov)
    fprintf (fp, "rhov = %g\n", tp->rhov (ts));
  if (tp->muv)
    fprintf (fp, "muv = %g\n", tp->muv (ts));
  if (tp->lambdav)
    fprintf (fp, "lambdav = %g\n", tp->lambdav (ts));
  if (tp->cpv)
    fprintf (fp, "cpv = %g\n", tp->cpv (ts));
  if (tp->pvap) {
    for (int jj=0; jj<NS; jj++)
      fprintf (fp, "pvap[%d] = %g\n", jj, tp->pvap (ts, jj));
  }
  if (tp->dhev) {
    for (int jj=0; jj<NS; jj++)
      fprintf (fp, "dhev[%d] = %g\n", jj, tp->dhev (ts, jj));
  }
  if (tp->diff) {
    for (int jj=0; jj<NS; jj++)
      fprintf (fp, "diff[%d] = %g\n", jj, tp->diff (ts, jj));
  }
  if (tp->cps) {
    for (int jj=0; jj<NS; jj++)
      fprintf (fp, "cps[%d] = %g\n", jj, tp->cps (ts, jj));
  }
  if (tp->sigmas) {
    for (int jj=0; jj<NS; jj++)
      fprintf (fp, "sigmas[%d] = %g\n", jj, tp->sigmas (ts, jj));
  }
}

/**
### *gasprop_thermal_expansion()*: Thermal expansion coefficient of an ideal gas
*/

double gasprop_thermal_expansion (ThermoState * ts) {
  return ts->T > 0. ? 1./ts->T : 0.;
}

/**
### *liqprop_thermal_expansion()*: Thermal expansion coefficient of a liquid
*/

double liqprop_thermal_expansion (ThermoProps * tp, ThermoState * ts) {
  double epsT = 1.e-3;
  double Ttop = ts->T + epsT, Tbot = ts->T - epsT;
  ThermoState tstop, tsbot;
  tstop.T = Ttop, tstop.P = ts->P, tstop.x = ts->x;
  tsbot.T = Tbot, tsbot.P = ts->P, tsbot.x = ts->x;
  double rhotop = tp->rhov (&tstop), rhobot = tp->rhov (&tsbot);
  double rhoval = tp->rhov (ts);
  return (rhoval > 0.) ? -1./rhoval*(rhotop - rhobot)/(2.*epsT) : 0.;
}

