#include "fractions.h"
#include "phase.h"
#include "vof.h"

#define rho(f) rho2v[]
scalar rho2v[], mu2v[];

scalar f[], * interfaces = {f};

face vector alphav[], muv[];
scalar rhov[];

int NS = 0;
double rhoval = 1., muval = 0., lambdaval = 1., cpval = 1., Dmixval = 1.;
double T0 = 300., YG0 = 1., Pref = 101325., MWval = 1.;

Phase * gas;
char ** gas_species = NULL;

struct PhaseChangeModel {
  bool isothermal;
  bool isomassfrac;
  bool fick_corrected;
  bool molar_diffusion;
  bool mass_diffusion_enthalpy;
  bool divergence;
  bool chemistry;
  bool normalize;
  bool no_advection_div;
} pcm = {
  false,    // isothermal
  false,    // isomassfrac
  true,     // fick_corrected
  true,     // molar_diffusion
  true,     // mass_diffusion_enthalpy
  true,     // divergence
  true,     // chemistry
  true,     // normalize
  true,     // no_advection_div
};

static scalar * f_tracers = NULL;

event defaults (i = 0) {
  gas = new_phase ("", NS, true, gas_species);

  gas->isothermal = pcm.isothermal;
  gas->isomassfrac = pcm.isomassfrac;

  double * xg = malloc (gas->n*sizeof (double));
  double * Dg = malloc (gas->n*sizeof (double));
  double * dhevsg = malloc (gas->n*sizeof (double));
  double * cpsg = malloc (gas->n*sizeof (double));
  double * MWg = malloc (gas->n*sizeof (double));

  foreach_species_in (gas)
    xg[i] = 1.;
  correctfrac (xg, gas->n);

  if (NS == 1)
    xg[0] = 1.;
  else if (NS == 2)
    xg[0] = YG0, xg[1] = 1. - YG0;

  ThermoState tsg;
  tsg.T = T0, tsg.P = Pref, tsg.x = xg;

  phase_set_thermo_state (gas, &tsg);

  foreach_species_in (gas) {
    Dg[i] = Dmixval;
    cpsg[i] = cpval;
    MWg[i] = MWval;
  }

  phase_set_properties (gas,
      rho = rhoval, mu = muval,
      lambda = lambdaval, cp = cpval,
      D = Dg, cps = cpsg, MWs = MWg);

  free (xg);
  free (Dg);
  free (dhevsg);
  free (cpsg);
  free (MWg);

  phase_set_tracers (gas);

  f_tracers = f.tracers;
  f.tracers = list_concat (f.tracers, gas->tracers);

  alpha = alphav;
  rho = rhov;
  mu = muv;

  rho2v.nodump = mu2v.nodump = true;

  no_advection_div = pcm.no_advection_div;
}

event cleanup (t = end) {
  free (f.tracers), f.tracers = f_tracers;
  delete_phase (gas);
}

event reset_sources (i++) {
  phase_reset_sources (gas);

  foreach() {
    for (scalar drhodt in drhodtlist)
      drhodt[] = 0.;
  }
}

event end_init (i = 0);

event reset_sources (i++);

#if VARIABLE_PROPERTIES
event phase_properties (i++) {
  phase_update_mw_moles (gas, f, tol = P_ERR, extend = true);
  phase_update_properties (gas, &tp2, f, P_ERR);
}
#endif

event chemistry (i++) {
#if CHEMISTRY
  if (pcm.chemistry) {
    ode_function batch = batch_isothermal_constantpressure;
    unsigned int NEQ = gas->n;
    if (!gas->isothermal) {
      batch = batch_nonisothermal_constantpressure;
      NEQ++;
    }
# if BINNING
    scalar Y = gas->YList[index_species ("H2")];
    scalar T = gas->T;
    phase_chemistry_binning (gas, dt, batch, NEQ, {T,Y}, (double[]){1e-2,1e-3},
        verbose = true, f = f, tol = 1-F_ERR);
# else
    phase_chemistry_direct (gas, dt, batch, NEQ, f, tol = 1-F_ERR);
# endif
  }
#endif

  if (pcm.mass_diffusion_enthalpy)
    phase_add_heat_species_diffusion (gas, f, pcm.molar_diffusion);
}

#if VARIABLE_PROPERTIES
event divergence (i++) {
  if (pcm.divergence) {
    phase_update_divergence (gas, f, pcm.fick_corrected, pcm.molar_diffusion);

    scalar divu2 = gas->divu;
    foreach()
      drhodt[] = divu2[];
  }
}
#endif

event tracer_advection (i++);

event tracer_diffusion (i++) {
  phase_update_mw_moles (gas, f, tol = P_ERR);
  phase_diffusion_velocity (gas, f, pcm.fick_corrected, pcm.molar_diffusion);
  phase_diffusion (gas, f, varcoeff = true);

  if (pcm.normalize)
    phase_normalize_fractions (gas);
}

event properties (i++) {
  foreach_face() {
    alphav.x[] = fm.x[]/face_value (rho2v, 0);

    face vector muv = mu;
    muv.x[] = fm.x[]*face_value (mu2v, 0);
  }

  foreach()
    rhov[] = cm[]*rho2v[];
}

event properties (i++) {
  scalar rhog = gas->rho;
  scalar mug = gas->mu;
  foreach() {
    rho2v[] = rhog[];
    mu2v[] = mug[];
  }
}

double CFL_MAX = 0.5;

event stability (i++) {
  if (CFL > CFL_MAX)
    CFL = CFL_MAX;
}

