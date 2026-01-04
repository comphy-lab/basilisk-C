/**
# Phase Model

We gather information about a multicomponent, non-isothermal, variable material
properties phase into a structure. This organization facilitates the solution of
advection-diffusion-reaction equations for each scalar field of the phase,
improving re-usability, and limiting code duplication. */

#include "variable-properties.h"
#include "common-evaporation.h"
#include "diffusion.h"

#define UNDEFINED -1
#ifndef F_ERR
# define F_ERR 1e-10
#endif

/**
## Phase structure

The following structure contains the attributes of the phase
*/

typedef struct {
  // Public attributes
  char * name;
  char ** species;
  scalar T;
  scalar P;
  scalar STimp;
  scalar STexp;
  scalar * YList;
  scalar * XList;
  scalar * SYimpList;
  scalar * SYexpList;
  size_t n;
  bool isothermal;
  bool isomassfrac;
  bool inverse;
  bool dump_all;
  scalar * tracers;
  // Material properties
  scalar rho;
  scalar mu;
  scalar MW;
  scalar lambda;
  scalar cp;
  scalar dhev;
  scalar * DList;
  scalar * cpList;
  scalar * dhevList;
  double * MWs;
  ThermoState * ts0;
  // Flow divergence
  scalar divu;
  scalar betaT;
  scalar * betaYList;
  scalar DTDt;
  scalar * DYDtList;
} Phase;

/**
## Iterators

The following iterators can be used to expose the scalar fields of the phase
and to loop over the species and relative properties. */

macro foreach_scalar_in (Phase * phase) {
  {
    scalar T = phase->T; NOT_UNUSED (T);
    scalar P = phase->P; NOT_UNUSED (P);
    scalar rho = phase->rho; NOT_UNUSED (rho);
    scalar mu = phase->mu; NOT_UNUSED (mu);
    scalar MW = phase->MW; NOT_UNUSED (MW);
    scalar lambda = phase->lambda; NOT_UNUSED (lambda);
    scalar cp = phase->cp; NOT_UNUSED (cp);
    scalar dhev = phase->dhev; NOT_UNUSED (dhev);
    scalar STexp = phase->STexp; NOT_UNUSED (STexp);
    scalar STimp = phase->STimp; NOT_UNUSED (STimp);
    scalar divu = phase->divu; NOT_UNUSED (divu);
    scalar betaT = phase->betaT; NOT_UNUSED (betaT);
    scalar DTDt = phase->DTDt; NOT_UNUSED (DTDt);
    {...}
  }
}

macro foreach_species_in (Phase * phase) {
  for (size_t i = 0; i < phase->n; i++) {
    scalar Y = {-1}; NOT_UNUSED (Y);
    scalar X = {-1}; NOT_UNUSED (X);
    scalar D = {-1}; NOT_UNUSED (D);
    scalar cps = {-1}; NOT_UNUSED (cps);
    scalar dhevs = {-1}; NOT_UNUSED (dhevs);
    scalar SYexp = {-1}; NOT_UNUSED (SYexp);
    scalar SYimp = {-1}; NOT_UNUSED (SYimp);
    scalar betaY = {-1}; NOT_UNUSED (betaY);
    scalar DYDt = {-1}; NOT_UNUSED (DYDt);
    if (phase->YList) Y = phase->YList[i];
    if (phase->XList) X = phase->XList[i];
    if (phase->DList) D = phase->DList[i];
    if (phase->cpList) cps = phase->cpList[i];
    if (phase->dhevList) dhevs = phase->dhevList[i];
    if (phase->SYexpList) SYexp = phase->SYexpList[i];
    if (phase->SYimpList) SYimp = phase->SYimpList[i];
    if (phase->betaYList) betaY = phase->betaYList[i];
    if (phase->DYDtList) DYDt = phase->DYDtList[i];
    {...}
  }
}

/**
## Constructors

We define macros that facilitates the creation of new scalar or vector fields
*/

#define new_field_scalar(Y, phase, no_dump)                                   \
  {                                                                           \
    scalar Y = new scalar;                                                    \
    phase->Y = Y;                                                             \
    Y.inverse = phase->inverse;                                               \
    char name[80];                                                            \
    sprintf (name, "%s%s", #Y, phase->name);                                  \
    free (Y.name);                                                            \
    Y.name = strdup (name);                                                   \
    Y.nodump = no_dump;                                                       \
  }

#define new_field_type(type, Y, phase, no_dump)                               \
  new_field_##type(Y, phase, no_dump);

#define scalar_name(name, Y, phase, i)                                        \
  sprintf (name, "%s%s%s", #Y, phase->name, phase->species[i]);               \

#define vector_name(name, Y, phase, i, ext)                                   \
  sprintf (name, "%s%s%s%s", #Y, phase->name, phase->species[i], ext.x)       \

#define new_field_scalar_name(Y, phase, no_dump)                              \
  scalar Y = new scalar;                                                      \
  Y.inverse = phase->inverse;                                                 \
  char name[80];                                                              \
  scalar_name (name, Y, phase, i);                                            \
  free (Y.name);                                                              \
  Y.name = strdup (name);                                                     \
  Y.nodump = no_dump;

#define new_list_scalar_name(Y, phase, list, no_dump)                         \
  for (size_t i = 0; i < phase->n; i++) {                                     \
    new_field_scalar_name (Y, phase, no_dump);                                \
    list = list_add (list, Y);                                                \
  }

#define new_field_vector_name(Y, phase, no_dump)                              \
  vector Y = new vector;                                                      \
  char name[80];                                                              \
  foreach_dimension() {                                                       \
    Y.x.inverse = phase->inverse;                                             \
    vector_name (name, Y, phase, i, ext);                                     \
    free (Y.x.name);                                                          \
    Y.x.name = strdup (name);                                                 \
    Y.x.nodump = no_dump;                                                     \
  }                                                                           \

#define new_list_vector_name(Y, phase, list, no_dump)                         \
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};                    \
  for (size_t i = 0; i < phase->n; i++) {                                     \
    new_field_vector_name(Y, phase, no_dump);                                 \
    list = vectors_add (list, Y);                                             \
  }

#define new_list_type_name(type, Y, phase, list, no_dump)                     \
  list = NULL;                                                                \
  new_list_##type##_name(Y, phase, list, no_dump)

/**
### *phase_species_names()*: populate species names with user provided names */

void phase_species_names (Phase * phase, char ** names = NULL) {
  if (names) {
    phase->species = (char **)malloc (phase->n*sizeof (char *));
    for (size_t i = 0; i < phase->n; i++)
      phase->species[i] = strdup (names[i]);
  }
  else {
    phase->species = (char **)malloc (phase->n*sizeof (char *));
    for (size_t i = 0; i < phase->n; i++) {
      char name[80];
      sprintf (name, "%zu", i+1);
      phase->species[i] = strdup (name);
    }
  }
}

/**
### *new_phase_empty()*: create an empty phase, i.e. without fields */

Phase * new_phase_empty (char * name = "", bool inverse = false) {
  Phase * phase = (Phase *)malloc (sizeof (Phase));
  phase->name = strdup (name);
  phase->n = 0;
  phase->YList = NULL;
  phase->XList = NULL;
  phase->MWs = NULL;
  phase->ts0 = NULL;
  phase->isothermal = true;
  phase->isomassfrac = true;
  phase->inverse = inverse;
  phase->tracers = NULL;
  phase->dump_all = false;

  /**
  All fields are undefined by default, while all lists are NULL. This is
  necessary to avoid segmentation from the foreach_species loops. */

  phase->T.i = -1;
  phase->P.i = -1;
  phase->rho.i = -1;
  phase->mu.i = -1;
  phase->MW.i = -1;
  phase->lambda.i = -1;
  phase->cp.i = -1;
  phase->dhev.i = -1;
  phase->STexp.i = -1;
  phase->STimp.i = -1;
  phase->divu.i = -1;
  phase->betaT.i = -1;
  phase->DTDt.i = -1;

  phase->YList = NULL;
  phase->XList = NULL;
  phase->DList = NULL;
  phase->cpList = NULL;
  phase->dhevList = NULL;
  phase->SYexpList = NULL;
  phase->SYimpList = NULL;
  phase->betaYList = NULL;
  phase->DYDtList = NULL;

  return phase;
}

/**
### *new_phase_minimal()*: create a minimal phase, i.e. without properties */

Phase * new_phase_minimal (char * name = "", size_t ns = 0,
    bool inverse = false, char ** species = NULL)
{
  Phase * phase = new_phase_empty (name, inverse);
  phase->n = ns;

  // Default set of names
  phase_species_names (phase, species);

  // Which fields must be dumped
  bool nodump = !phase->dump_all;

  // Create minimal set of scalar fields
  new_list_type_name (scalar, Y, phase, phase->YList, false);
  new_list_type_name (scalar, X, phase, phase->XList, nodump);
  new_field_type (scalar, T, phase, false);

  // Create minimal set of material properties
  new_field_type (scalar, MW, phase, false);

  foreach() {
    foreach_scalar_in (phase) {
      T[] = 0.;
      MW[] = 0.;
      foreach_species_in (phase) {
        Y[] = 0.;
        X[] = 0.;
      }
    }
  }

  // Create species molecular weights and set it to 1
  phase->MWs = (double *)malloc (phase->n*sizeof (double));
  foreach_species_in (phase)
    phase->MWs[i] = 1.;

  return phase;
}

/**
### *new_phase()*: create a field with all fields of of the structure */

Phase * new_phase (char * name = "", size_t ns = 0, bool inverse = false,
    char ** species = NULL)
{
  Phase * phase = new_phase_empty (name, inverse);
  phase->n = ns;
  phase->isothermal = false;
  phase->isomassfrac = false;

  // Default set of names
  phase_species_names (phase, species);

  // Which fields must be dumped
  bool nodump = !phase->dump_all;

  // Create scalar fields
  new_field_type (scalar, T, phase, false);
  new_list_type_name (scalar, Y, phase, phase->YList, false);
  new_list_type_name (scalar, X, phase, phase->XList, nodump);

  // Create material properties
  new_field_type (scalar, P, phase, false);
  new_field_type (scalar, rho, phase, false);
  new_field_type (scalar, mu, phase, false);
  new_field_type (scalar, MW, phase, false);
  new_list_type_name (scalar, D, phase, phase->DList, false);
  new_list_type_name (scalar, cp, phase, phase->cpList, false);
  new_list_type_name (scalar, dhev, phase, phase->dhevList, false);
  new_field_type (scalar, lambda, phase, false);
  new_field_type (scalar, cp, phase, false);
  new_field_type (scalar, dhev, phase, false);
  new_field_type (scalar, divu, phase, nodump);
  new_field_type (scalar, betaT, phase, nodump);
  new_field_type (scalar, DTDt, phase, nodump);

  // Create source terms
  new_field_type (scalar, STimp, phase, nodump);
  new_field_type (scalar, STexp, phase, nodump);
  new_list_type_name (scalar, SYimp, phase, phase->SYimpList, nodump);
  new_list_type_name (scalar, SYexp, phase, phase->SYexpList, nodump);
  new_list_type_name (scalar, betaY, phase, phase->betaYList, nodump);
  new_list_type_name (scalar, DYDt, phase, phase->DYDtList, nodump);

  foreach()
    foreach_scalar_in (phase) {
      T[] = 0.;
      P[] = 0.;
      rho[] = 0.;
      mu[] = 0.;
      MW[] = 1.;
      lambda[] = 0.;
      cp[] = 0.;
      dhev[] = 0.;
      STimp[] = 0.;
      STexp[] = 0.;
      divu[] = 0.;
      betaT[] = 0.;
      DTDt[] = 0.;
      foreach_species_in (phase) {
        Y[] = 0.;
        X[] = 0.;
        D[] = 0.;
        cps[] = 0.;
        dhevs[] = 0.;
        SYimp[] = 0.;
        SYexp[] = 0.;
        betaY[] = 0.;
        DYDt[] = 0.;
      }
    }

  // Create species molecular weights and set it to 1
  phase->MWs = (double *)malloc (phase->n*sizeof (double));
  foreach_species_in (phase)
    phase->MWs[i] = 1.;

  // Create the initial thermo state
  phase->ts0 = new_thermo_state (phase->n);

#if TREE
  scalar rhov = phase->rho;
  rhov.restriction = density_restriction;
  rhov.refine = rhov.prolongation = density_refine;
  rhov.dirty = true;
#endif

  return phase;
}

/**
### *delete_phase()*: free phases from memory */

void delete_phase (Phase * phase) {
  foreach_scalar_in (phase)
    delete ({T,P,STimp,STexp,rho,mu,MW,lambda,cp,dhev,divu,betaT,DTDt});

  if (phase->YList) delete (phase->YList), free (phase->YList);
  if (phase->XList) delete (phase->XList), free (phase->XList);
  if (phase->SYimpList) delete (phase->SYimpList), free (phase->SYimpList);
  if (phase->SYexpList) delete (phase->SYexpList), free (phase->SYexpList);
  if (phase->DList) delete (phase->DList), free (phase->DList);
  if (phase->cpList) delete (phase->cpList), free (phase->cpList);
  if (phase->dhevList) delete (phase->dhevList), free (phase->dhevList);
  if (phase->betaYList) delete (phase->betaYList), free (phase->betaYList);
  if (phase->DYDtList) delete (phase->DYDtList), free (phase->DYDtList);

  if (phase->species)
    foreach_species_in (phase)
      free (phase->species[i]);

  if (phase->name) free (phase->name), phase->name = NULL;
  if (phase->species) free (phase->species), phase->species = NULL;
  if (phase->tracers) free (phase->tracers), phase->tracers = NULL;
  if (phase->ts0) free_thermo_state (phase->ts0);
  if (phase->MWs) free (phase->MWs);

  free (phase), phase = NULL;
}

/**
## Initial Conditions

The following functions can be used to set the initial conditions of the phase
*/

/**
### *phase_set_properties()*: set phase material property fields */

#define provided(x) (x != UNDEFINED)

#define provided_list(x) (x != NULL)

#define error_unprovided(message)                                             \
  {                                                                           \
    fprintf (ferr, "src/phase.h:%d: error: %s not provided\n",                \
        LINENO, #message);                                                    \
    fflush(ferr);                                                             \
    abort();                                                                  \
  }

#define check_provided(x)                                                     \
  if (!provided (x))                                                          \
    error_unprovided (x)

#define check_provided_list(x)                                                \
  if (!provided_list (x))                                                     \
    error_unprovided (x)

void phase_set_properties (Phase * phase,
    double rho = UNDEFINED, double mu = UNDEFINED,
    double lambda = UNDEFINED, double cp = UNDEFINED,
    double dhev = UNDEFINED, double * dhevs = NULL,
    double * D = NULL, double * MWs = NULL,
    double * cps = NULL)
{
  if (provided_list (MWs))
    foreach_species_in (phase)
      phase->MWs[i] = MWs[i];

  foreach() {
    if (provided (rho)) {
      scalar phase_rho = phase->rho;
      phase_rho[] = rho;
    }
    if (provided (mu)) {
      scalar phase_mu = phase->mu;
      phase_mu[] = mu;
    }
    if (provided (dhev)) {
      scalar phase_dhev = phase->dhev;
      phase_dhev[] = dhev;
    }

    if (!phase->isomassfrac) {
      for (size_t i = 0; i < phase->n; i++) {
        if (provided_list (D)) {
          scalar phase_D = phase->DList[i];
          phase_D[] = D[i];
        }
        if (provided_list (dhevs)) {
          scalar phase_dhev = phase->dhevList[i];
          phase_dhev[] = dhevs[i];
        }
        if (provided_list (cps)) {
          scalar phase_cp = phase->cpList[i];
          phase_cp[] = cps[i];
        }
      }
    }

    if (!phase->isothermal) {
      if (provided (lambda)) {
        scalar phase_lambda = phase->lambda;
        phase_lambda[] = lambda;
      }
      if (provided (cp)) {
        scalar phase_cp = phase->cp;
        phase_cp[] = cp;
      }
    }

  }
}

/**
### *phase_is_uniform()*: check if species and temperature are constant in space */

bool phase_is_uniform (Phase * phase) {
  bool uniform = true;
  foreach_scalar_in (phase) {
    double Tmin = statsf (T).min;
    double Tmax = statsf (T).max;
    uniform &= (Tmin == Tmax);
    foreach_species_in (phase) {
      double Ymin = statsf (Y).min;
      double Ymax = statsf (Y).max;
      uniform &= (Ymin == Ymax);
    }
  }
  return uniform;
}

/**
### *phase_species_index()*: index of a species given its name */

size_t phase_species_index (Phase * phase, char * species) {
  if (phase->species) {
    foreach_species_in (phase) {
      if (strcmp (phase->species[i], species) == 0)
        return i;
    }
    fprintf (ferr, "src/phase.h:%d: error: species %s not found\n",
        LINENO, species), fflush (ferr);
    abort();
  }
  else {
    fprintf (ferr, "src/phase.h:%d: error: species names not provided\n",
        LINENO), fflush (ferr);
    abort();
  }
}

/**
### *phase_set_thermo_state()*: set the initial thermodynamic state */

void phase_set_thermo_state (Phase * phase, const ThermoState * ts,
    bool force = false)
{
  copy_thermo_state (phase->ts0, ts, phase->n);
  if (phase_is_uniform (phase) || force) {
    foreach() {
      foreach_scalar_in (phase) {
        T[] = ts->T;
        P[] = ts->P;
        foreach_species_in (phase)
          if (ts->x) Y[] = ts->x[i];
      }
    }
  }
}

/**
### *phase_set_composition_from_string()*: set species mass fractions

Usage: phase_set_composition_from_string (phase, "NC7H16 0.2 N2 0.8"); */

void phase_set_composition_from_string (Phase * phase, char * s,
    char * sep = " ", bool force = false)
{
  ThermoState * ts = new_thermo_state (phase->n);
  ts->T = phase->ts0->T;
  ts->P = phase->ts0->P;
  foreach_species_in (phase)
    ts->x[i] = 0.;

  char * input = strdup (s);
  char * token = strtok (input, sep);
  unsigned int count = 0;
  while (token != NULL) {
    count++;
    char species[80];
    double val = 0.;
    if (count%2 != 0) {
      strcpy (species, token);
      val = 0;
    }
    else {
      val = atof (token);
      size_t index = phase_species_index (phase, species);
      ts->x[index] = val;
    }
    token = strtok (NULL, sep);
  }
  phase_set_thermo_state (phase, ts, force);

  free_thermo_state (ts);
  free (input);
}


/**
## Solvers and Utilities

Functions of the phase, which can be useful in different contexts. */

/**
### *phase_set_tracers()*: populate the tracers list for advection

The list `phase->tracers` can then be used to perform the advection using the
methods implemented in [tracer.h](src/tracer.h) or in [vof.h](src/vof.h). */

void phase_set_tracers (Phase * phase) {
  if (!phase->isomassfrac)
    phase->tracers = list_concat (phase->tracers, phase->YList);
  if (!phase->isothermal)
    phase->tracers = list_add (phase->tracers, phase->T);
}

/**
### *phase_set_gradient()*: set a flux limiter to each phase tracer */

void phase_set_gradient (Phase * phase,
    double (* gradient) (double, double, double))
{
  for (scalar s in phase->tracers)
    s.gradient = gradient;
}

/**
### *phase_reset_sources()*: set to zero source term fields */

void phase_reset_sources (Phase * phase) {
  foreach() {
    foreach_scalar_in (phase) {
      STexp[] = 0.;
      STimp[] = 0.;
      divu[] = 0.;
      DTDt[] = 0.;
      foreach_species_in (phase) {
        DYDt[] = 0.;
        SYexp[] = 0.;
        SYimp[] = 0.;
      }
    }
  }
}

/**
### *phase_diffusion()*: resolve diffusion of species and temperature */

void phase_diffusion (Phase * phase, (const) scalar f = unity,
    bool varcoeff = false)
{
  scalar ff[];
  foreach()
    ff[] = phase->inverse ? 1. - f[] : f[];

  face vector fs[];
  face_fraction (ff, fs); // fixme: can't use f in this function

  foreach_scalar_in (phase) {
    if (!phase->isothermal) {
      scalar thetaT[];
      foreach()
        thetaT[] = varcoeff ? cm[]*max (ff[]*rho[]*cp[], F_ERR)
                            : cm[]*max (ff[], F_ERR);

      face vector lambdaf[];
      foreach_face() {
        lambdaf.x[] = varcoeff ? face_value (lambda, 0)
          : face_value (lambda, 0) /
            (face_value (rho, 0)*face_value (cp, 0) + 1e-16);
        lambdaf.x[] *= fs.x[]*fm.x[];
      }

      if (!varcoeff)
        foreach() {
          STexp[] = (rho[]*cp[] > 0.) ? STexp[]/(rho[]*cp[]) : STexp[];
          STimp[] = (rho[]*cp[] > 0.) ? STimp[]/(rho[]*cp[]) : STimp[];
        }

      diffusion (T, dt, D=lambdaf, r=STexp, beta=STimp, theta=thetaT);
    }

    if (!phase->isomassfrac) {
      foreach_species_in (phase) {
        scalar thetaY[];
        foreach()
          thetaY[] = varcoeff ? cm[]*max (ff[]*rho[], F_ERR)
                              : cm[]*max (ff[], F_ERR);

        face vector Df[];
        foreach_face() {
          Df.x[] = varcoeff ? face_value (rho, 0)*face_value (D, 0)
                            : face_value (D, 0);
          Df.x[] *= fs.x[]*fm.x[];
        }

        if (!varcoeff)
          foreach() {
            SYexp[] = (rho[] > 0.) ? SYexp[]/rho[] : SYexp[];
            SYimp[] = (rho[] > 0.) ? SYimp[]/rho[] : SYimp[];
          }

        diffusion (Y, dt, D=Df, r=SYexp, beta=SYimp, theta=thetaY);
      }

    }
  }

  phase_reset_sources (phase);  // fixme: maybe I still need those terms
}

/**
### *phase_tracers_to_scalars()*: multiply each scalar by a volume fraction */

void phase_tracers_to_scalars (Phase * phase, scalar f, double tol = 1e-10) {
  foreach() {
    double ff = phase->inverse ? 1. - f[] : f[];
    foreach_scalar_in (phase) {
      if (!phase->isothermal)
        T[] = (ff > tol) ? T[]/ff : 0.;
      foreach_species_in (phase)
        if (!phase->isomassfrac)
          Y[] = (ff > tol) ? Y[]/ff : 0.;
    }
  }
}

/**
### *phase_scalars_to_tracers()*: divide each scalar by a volume fraction */

void phase_scalars_to_tracers (Phase * phase, scalar f) {
  foreach() {
    double ff = phase->inverse ? 1. - f[] : f[];
    foreach_scalar_in (phase) {
      if (!phase->isothermal)
        T[] *= phase->isothermal ? 1. : ff;
      foreach_species_in (phase)
        if (!phase->isomassfrac)
          Y[] *= ff;
    }
  }
}

/**
### *phase_normalize_fractions()*: normalize mass fractions to 1

When we resolve the diffusion step with variable diffusivity, we need to account
for explicit corrections which enforce mass conservation. These corrections are
necessary but they do not work perfectly, i.e. the sum of the mass fractions
will not be exactly 1. This brutal normalization helps respecting this
constraint. */

void phase_normalize_fractions (Phase * phase) {
  double * ymass = (double *)malloc (phase->n*sizeof (double));
  foreach(serial) {
    foreach_species_in (phase)
      ymass[i] = Y[];
    correctfrac (ymass, phase->n);
    foreach_species_in (phase)
      Y[] = ymass[i];
  }
  free (ymass);
}

/**
### *phase_add_heat_species_diffusion()*: heat contribution from mass diffusion

This term comes from the energy equation for a multicomponent mixture, solved
for the temperature field. */

void phase_add_heat_species_diffusion (Phase * phase, (const) scalar f = unity,
    bool molar_diffusion = false, double tol = 1e-10)
{
  if (!phase->isothermal) {
    foreach_scalar_in (phase) {
      foreach() {
        double mde = 0.;
        coord gT = {0.,0.,0.};
        coord gY = {0.,0.,0.};
        coord gYsum = {0.,0.,0.};

        foreach_dimension()
          gT.x = (T[1] - T[-1])/(2.*Delta);

        foreach_dimension() {
          foreach_species_in (phase) {
            if (molar_diffusion)
              gYsum.x -= (MW[] > 0.) ?
                rho[]*D[]*phase->MWs[i]/MW[]*(X[1] - X[-1])/(2.*Delta) : 0.;
            else
              gYsum.x -= rho[]*D[]*(Y[1] - Y[-1])/(2.*Delta);
          }

          foreach_species_in (phase) {
            if (molar_diffusion)
              gY.x = (MW[] > 0.) ?
                -rho[]*D[]*phase->MWs[i]/MW[]*(X[1] - X[-1])/(2.*Delta) : 0.;
            else
              gY.x = -rho[]*D[]*(Y[1] - Y[-1])/(2.*Delta);
            mde += cps[]*(gY.x - Y[]*gYsum.x)*gT.x;
          }
        }
        double ff = phase->inverse ? 1. - f[] : f[];
        STexp[] -= mde*cm[]*(ff == 1);  // fixme: use interfacial gradients if
                                        // you want to apply it to interfacial
                                        // cells
      }
    }
  }
}

/**
## Thermodynamic, transport, and kinetics

The following functions allow including the effect of variable material
properties and the chemistry integration. */

/**
### *phase_update_mw_moles()*: compute mole fractions and molecular weights

We resolve for the mass fractions and update the mole fractions at needs, as
they are used for the variable material properties calculation. The mixture
molecular weight is also computed during this step. */

void phase_update_mw_moles (Phase * phase, (const) scalar f = unity,
    double tol = 1e-10, bool extend = false)
{
  double * ymass = (double *)malloc (phase->n*sizeof (double));
  double * xmole = (double *)malloc (phase->n*sizeof (double));
  foreach_scalar_in (phase) {
    foreach(serial) {
      //double ff = phase->inverse ? 1. - f[] : f[];
      //if (ff > tol) {
        foreach_species_in (phase)
          ymass[i] = (phase->n == 1.) ? 1. : Y[];
        correctfrac (ymass, phase->n);
        mass2molefrac (xmole, ymass, phase->MWs, phase->n);
        MW[] = mass2mw (ymass, phase->MWs, phase->n);
        foreach_species_in (phase)
          X[] = xmole[i];
      //}
    }

    if (extend) {
      foreach(serial) {
        double ff = phase->inverse ? 1. - f[] : f[];
        if (ff <= tol) {
          int counter = 0;
          double ext_MW = 0.;
          foreach_neighbor(1) {
            double ffnei = phase->inverse ? 1. - f[] : f[];
            if (ffnei > tol) {
              counter++;
              ext_MW += MW[];
            }
          }
          MW[] = (counter != 0.) ? ext_MW/counter : 0.;
        }
      }
    }
  }
  free (ymass);
  free (xmole);
}

/**
### *phase_update_properties()*: update material properties

Compute variable material property fields, according to the thermodynamic state
(temperature, pressure, and composition, and using the functions of the
thermodynamic properties structure. */

void phase_update_properties (Phase * phase, const ThermoProps * tp,
    (const) scalar f = unity, double tol = 1e-10)
{
  double * xmole = (double *)malloc (phase->n*sizeof (double));
  double * arrdiff = (double *)malloc (phase->n*sizeof (double));
  double * arrbetaY = (double *)malloc (phase->n*sizeof (double));
  double * arrdhev = (double *)malloc (phase->n*sizeof (double));
  double * arrcps  = (double *)malloc (phase->n*sizeof (double));
  foreach(serial) {
    double ff = phase->inverse ? 1. - f[] : f[];
    if (ff > tol) {
      ThermoState ts;
      foreach_scalar_in (phase) {
        foreach_species_in (phase)
          xmole[i] = X[];
        ts.T = T[];
        ts.P = P[];
        ts.x = (phase->n == 1) ? (double[]){1.} : xmole;

        if (tp->rhov) rho[] = tp->rhov (&ts);
        if (tp->muv) mu[] = tp->muv (&ts);
        if (tp->lambdav) lambda[] = tp->lambdav (&ts);
        if (tp->cpv) cp[] = tp->cpv (&ts);
         if (tp->betaT) betaT[] = tp->betaT (tp, &ts);

        if (tp->diff) tp->diff (&ts, arrdiff);
        if (tp->betaY) tp->betaY (tp, &ts, arrbetaY);
        if (tp->dhev) tp->dhev (&ts, arrdhev);
        if (tp->cps)  tp->cps  (&ts, arrcps);

        foreach_species_in (phase) {
          if (tp->diff) D[] = arrdiff[i];
          if (tp->betaY) betaY[] = arrbetaY[i];
          if (tp->dhev) dhevs[] = arrdhev[i];
          if (tp->cps) cps[] = arrcps[i];
        }
      }
    }
  }
  free (xmole);
  free (arrdiff);
  free (arrbetaY);
  free (arrdhev);
  free (arrcps);
}

/**
### *phase_extend_properties()*: extend properties across the interface

To avoid SIGFPE problems in dry cells we extend the variable property fields
around the interface (using a brutal averaging approach). */

#define increment_property(ext_s, s) \
  ext_s += (s.i > 0) ? s[] : 0;

#define assign_property(ext_s, s, counter) \
  if (s.i > 0) s[] = (counter > 0) ? ext_s/counter : 0.;

void phase_extend_properties (Phase * phase,
    (const) scalar f = unity, double tol = 1e-10)
{
  double * ext_diff = (double *)malloc (phase->n*sizeof (double));
  double * ext_cps = (double *)malloc (phase->n*sizeof (double));
  double * ext_dhevs = (double *)malloc (phase->n*sizeof (double));
  double * ext_betaY = (double *)malloc (phase->n*sizeof (double));

  foreach_scalar_in (phase) {
    foreach(serial) {
      double ff = phase->inverse ? 1. - f[] : f[];
      if (ff <= tol) {
        double ext_rho = 0.;
        double ext_mu = 0.;
        double ext_MW = 0.;
        double ext_lambda = 0.;
        double ext_cp = 0.;
        double ext_dhev = 0.;
        double ext_betaT = 0.;

        foreach_species_in (phase) {
          ext_diff[i] = 0.;
          ext_cps[i] = 0.;
          ext_dhevs[i] = 0.;
          ext_betaY[i] = 0.;
        }

        int counter = 0;
        foreach_neighbor(1) {
          double ffnei = phase->inverse ? 1. - f[] : f[];
          if (ffnei > tol) {
            counter++;
            increment_property (ext_rho, rho);
            increment_property (ext_mu, mu);
            increment_property (ext_MW, MW);
            increment_property (ext_lambda, lambda);
            increment_property (ext_cp, cp);
            increment_property (ext_dhev, dhev);
            increment_property (ext_betaT, betaT);

            foreach_species_in (phase) {
              increment_property (ext_diff[i], D);
              increment_property (ext_cps[i], cps);
              increment_property (ext_dhevs[i], dhevs);
              increment_property (ext_betaY[i], betaY);
            }
          }
        }
        assign_property (ext_rho, rho, counter);
        assign_property (ext_mu, mu, counter);
        assign_property (ext_MW, MW, counter);
        assign_property (ext_lambda, lambda, counter);
        assign_property (ext_cp, cp, counter);
        assign_property (ext_dhev, dhev, counter);
        assign_property (ext_betaT, betaT, counter);

        foreach_species_in (phase) {
          assign_property (ext_diff[i], D, counter);
          assign_property (ext_cps[i], cps, counter);
          assign_property (ext_dhevs[i], dhevs, counter);
          assign_property (ext_betaY[i], betaY, counter);
        }
      }
    }
  }
  free (ext_diff);
  free (ext_cps);
  free (ext_dhevs);
  free (ext_betaY);
}

/**
### *phase_update_divergence()*: compute velocity divergence for this phase */

void phase_update_divergence (Phase * phase,
    (const) scalar f = unity,
    bool fick_corrected = true, bool molar_diffusion = true)
{
  scalar ff[];
  foreach()
    ff[] = phase->inverse ? 1. - f[] : f[];

  face vector fs[];
  face_fraction (ff, fs); // fixme: can't use f in this function

  foreach_scalar_in (phase) {

    /**
    We calculate the Lagrangian derivative of the temperature field. */

    face vector lambdagrad[];
    foreach_face()
      lambdagrad.x[] = face_value (lambda, 0)*face_gradient_x (T, 0) *
        fm.x[]*fs.x[];

    foreach() {
      foreach_dimension()
        DTDt[] += (lambdagrad.x[1] - lambdagrad.x[])/Delta;
      DTDt[] += STexp[] + STimp[]*T[];
    }

    /**
    We calculate the Lagrangian derivative of the chemical species mass fractions.
    */

    foreach_species_in (phase) {
      face vector rhoDmixY[];
      foreach_face() {
        double rhof = face_value (rho, 0);
        double Df = face_value (D, 0);
        rhoDmixY.x[] = rhof*Df*face_gradient_x (Y, 0)*fm.x[]*fs.x[];
      }

      foreach() {
        foreach_dimension()
          DYDt[] += (rhoDmixY.x[1] - rhoDmixY.x[])/Delta;
        DYDt[] += (SYexp[] + SYimp[]*Y[]); // fixme: I was using the YGInt
      }
    }

    /**
    We add diffusion correction contribution to the chemical species mass
    fraction derivatives. */
#if 1
    face vector phictot[];
    foreach_face() {
      phictot.x[] = 0.;
      if (fick_corrected) {
        foreach_species_in (phase) {
          double rhof = face_value (rho, 0);
          double Df = face_value (D, 0);
          if (molar_diffusion) {
            double MWf = face_value (MW, 0);
            phictot.x[] += (MWf > 0.) ?
              rhof*Df*phase->MWs[i]/MWf*face_gradient_x (X, 0)*fm.x[]*fs.x[] :
              0.;
          }
          else {
            phictot.x[] += rhof*Df*face_gradient_x (Y, 0)*fm.x[]*fs.x[];
          }
        }
      }
    }

    foreach_species_in (phase) {
      face vector phic[];
      foreach_face() {
        phic.x[] = phictot.x[];
        if (molar_diffusion) {
          double rhof = face_value (rho, 0);
          double Df = face_value (D, 0);
          double MWf = face_value (MW, 0);

          phic.x[] -= (MWf > 0.) ?
            rhof*Df/MWf*face_gradient_x (MW, 0)*fm.x[]*fs.x[] : 0.;
        }
        phic.x[] *= face_value (Y, 0);
      }

      foreach()
        foreach_dimension()
          DYDt[] -= (phic.x[1] - phic.x[])/Delta;
    }
#endif
    /**
    We calculate the divergence of the phase, just in the region occupied by the
    phase, defined by the volume fraction `f`. */

    foreach() {
      divu[] = 0.;

      // Add temperature contribution
      divu[] += (rho[]*cp[] > 0.) ? betaT[]/(rho[]*cp[])*DTDt[] : 0.;

      // Add chemical species contribution
      double divuspecies = 0.;
      foreach_species_in (phase)
        divuspecies += (rho[] > 0.) ? betaY[]/rho[]*DYDt[] : 0.;
      divu[] += divuspecies;

      // Volume-averaged divergence
      divu[] *= ff[];

      // Adjust sign for internal convention
      divu[] *= -1;
    }

  }
}

/**
### *phase_update_divergence_density()*: compute velocity divergence

Same as the previous function, but trying to discretize the density Lagrangian
derivative directly. This function has not been tested, and in principle it
should use the face velocity rather than the centered one. */

#if 0
void phase_update_divergence_density (Phase * phase, vector u,
    (const) scalar f = unity)
{
  foreach_scalar_in (phase) {
    vector grho[];
    gradients ({rho}, {grho});

    foreach() {
      divu[] = (rho[] - rho0[])/dt;

      foreach_dimension()
        divu[] += u.x[]*grho.x[];

      double ff = phase->inverse ? 1. - f[] : f[];
      divu[] *= (rho[] > 0.) ? cm[]/rho[]*ff : 0.;
    }
  }
}
#endif

/**
### *phase_diffusion_velocity()*: diffusion correction step

When using variable diffusivity (or different diffusivity values for each
species), Fick's law does not respect mass conservation. Therefore, we sum up
the error on the mass conservation, and we introduce it as an additional source
term of the equation, weighted on the mass fraction of each species
(`fick_corrected`). Additionally, the `molar_diffusion` correction accounts for
the fact that the diffusivity values may be computed on a mole-fraction-based
approach. */

void phase_diffusion_velocity (Phase * phase, (const) scalar f = unity,
    bool fick_corrected = true, bool molar_diffusion = true)
{
  if (!phase->isomassfrac) {
    scalar ff[];
    foreach()
      ff[] = phase->inverse ? 1. - f[] : f[];

    face vector fs[];
    face_fraction (ff, fs); // fixme: can't use f in this function

    foreach_scalar_in (phase) {
      face vector phictot[];
      foreach_face() {
        phictot.x[] = 0.;
        if (fick_corrected) {
          foreach_species_in (phase) {
            double rhof = face_value (rho, 0);
            double Df = face_value (D, 0);
            if (molar_diffusion) {
              double MWf = face_value (MW, 0);
              phictot.x[] += (MWf > 0.) ?
                rhof*Df*phase->MWs[i]/MWf*face_gradient_x (X, 0) : 0.;
            }
            else
              phictot.x[] += rhof*Df*face_gradient_x (Y, 0);
          }
        }
      }

      foreach_species_in (phase) {
        face vector phic[];
        foreach_face() {
          phic.x[] = phictot.x[];
          if (molar_diffusion) {
            double rhof = face_value (rho, 0);
            double Df = face_value (D, 0);
            double MWf = face_value (MW, 0);
            phic.x[] -= (MWf > 0.) ? rhof*Df/MWf*face_gradient_x (MW, 0) : 0.;
          }
          phic.x[] *= fm.x[]*fs.x[];
        }
        //double (* gradient_backup)(double, double, double) = Y.gradient;
        //Y.gradient = NULL;
        //face vector flux[];
        //tracer_fluxes (Y, phic, flux, dt, zeroc);
        //Y.gradient = gradient_backup;

        //foreach()
        //  foreach_dimension()
        //    SYexp[] +=(flux.x[] - flux.x[1])/Delta;

        face vector flux[];
        foreach_face()
          flux.x[] = face_value (Y, 0)*phic.x[];

#if !EMBED
        foreach()
          foreach_dimension()
            Y[] += (rho[] > 0.) ? dt/rho[]*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.;
#else
        update_tracer (f, phic, flux, dt);
#endif

#if 0
        foreach()
          foreach_dimension()
            SYexp[] += (flux.x[] - flux.x[1])/Delta;

        foreach()
          foreach_species_in (phase)
            DYDt[] += (flux.x[] - flux.x[1])/Delta;
#endif
      }
    }
  }
}

/**
### *phase_chemistry_direct()* direct integration of chemical reactions */

#if CHEMISTRY
void phase_chemistry_direct (Phase * phase, double dt,
    ode_function batch, unsigned int NEQ,
    (const) scalar f = unity, double tol = 1e-10)
{
  double * y0 = (double *)malloc (NEQ*sizeof (double));
  double * s0 = (double *)malloc (NEQ*sizeof (double));

  foreach_scalar_in (phase) {
    foreach(serial) {
      double ff = phase->inverse ? 1. - f[] : f[];
      if (ff > tol) {

        // Gather initial conditions
        foreach_species_in (phase)
          y0[i] = Y[];
        if (!phase->isothermal)
          y0[phase->n] = T[];

        // Additional data to be passed to the ODE function
        UserDataODE data;
        data.rho = rho[];
        data.cp = cp[];
        data.P = P[];
        data.T = T[];
        data.sources = s0;

        // Resolve the ODE system
        stiff_ode_solver (batch, NEQ, dt, y0, &data);

        // Recover the results of the ODE system
        foreach_species_in (phase) {
          Y[] = y0[i];
          DYDt[] += s0[i]*cm[];
        }
        if (!phase->isothermal) {
          T[] = y0[phase->n];
          DTDt[] += s0[phase->n]*cm[];
        }
      }
    }
  }
  free (y0);
  free (s0);
}

/**
### *phase_chemistry_binning()*: chemistry integration with binning */

# if BINNING
#include "binning.h"

scalar BINID[];

// fixme: to test
void phase_chemistry_binning (Phase * phase, double dt,
    ode_function batch, unsigned int NEQ,
    scalar * targets, double * eps, bool verbose = false,
    (const) scalar f = unity, double tol = 1e-10)
{
  double * s0 = (double *)malloc (NEQ*sizeof (double));

  scalar mask[];
  foreach() {
    double ff = phase->inverse ? 1. - f[] : f[];
    mask[] = (ff > tol) ? 1 : 0;
  }

  assert (NEQ == list_len (phase->tracers));
  scalar * fields = phase->tracers;

  // We store the old temperature and mass fractions for the calculation of the
  // divergence
  scalar * Y0List = list_clone (phase->YList);
  scalar T0[];
  foreach_scalar_in (phase) {
    foreach() {
      foreach_species_in (phase) {
        scalar Y0 = Y0List[i];
        Y0[] = Y[];
      }
      if (!phase->isothermal)
        T0[] = T[];
    }
  }

  // Split the domain in bins and return the table
  BinTable * table = binning (fields, targets, eps,
      phase->rho, phase->cp, mask = mask);

  // Bin-wise integration
  foreach_bin (table) {
    UserDataODE data;
    data.rho = bin_average (bin, phase->rho);
    data.cp = bin_average (bin, phase->cp);
    data.P = bin_average (bin, phase->P);
    data.T = bin_average (bin, phase->T);
    data.sources = s0;

    stiff_ode_solver (batch, NEQ, dt, bin->phi, &data);

    bin->rho = data.rho;
    bin->cp = data.cp;
  }
  binning_remap (table, fields, phase->rho, phase->cp);

  // Recover the source term for the divergence
  foreach_scalar_in (phase) {
    foreach() {
      foreach_species_in (phase) {
        scalar Y0 = Y0List[i];
        DYDt[] += (rho[] > 0) ? (Y[] - Y0[])/dt/rho[]*cm[] : 0.;
      }
      if (!phase->isothermal)
        DTDt[] += (rho[]*cp[] > 0) ? (T[] - T0[])/dt/rho[]/cp[]*cm[] : 0.;
    }
  }

  if (verbose) {
    bstats bs = binning_stats (table);
    fprintf (stdout, "bs->nactive = %zu\n", bs.nactive);
    fprintf (stdout, "bs->navg    = %zu\n", bs.navg);
    fprintf (stdout, "bs->nmax    = %zu\n", bs.nmax);
    fprintf (stdout, "bs->nmin    = %zu\n", bs.nmin);
    fprintf (stdout, "bs->nmask   = %zu\n", bs.nmask);
    fprintf (stdout, "\n");
  }
  binning_ids (table, BINID);

  binning_cleanup (table);
  delete (Y0List), free (Y0List);
  free (s0);
}
# endif   // BINNING
#endif    // CHEMISTRY

/**
## Phase Mass Balance

The following structure and function simplifies setting up mass balances for the
phase components. */

typedef struct {
  double * m;       // species mass in the domain
  double * m0;      // initial species mass in the domain
  double mtot;      // total mass in the domain
  double mtot0;     // initial mass in the domain
  double * mevap;   // species evaporated mass
  double mevaptot;  // total mass evaporated
  double * mf;      // species flux though the domain
  double mftot;     // total mass through the domain
  bool first;       // flag the first iteration
} PhaseMassBalance;

/**
### *new_phase_mass_balance()*: create a mass balance object */

PhaseMassBalance * new_phase_mass_balance (const Phase * phase) {
  PhaseMassBalance * pmb = malloc (sizeof (PhaseMassBalance));
  pmb->mtot = 0.;
  pmb->mtot0 = 0.;
  pmb->mftot = 0.;
  pmb->mevaptot = 0.;
  pmb->m = malloc (phase->n*sizeof (double));
  pmb->m0 = malloc (phase->n*sizeof (double));
  pmb->mevap = malloc (phase->n*sizeof (double));
  pmb->mf = malloc (phase->n*sizeof (double));
  foreach_species_in (phase) {
    pmb->m[i] = 0.;
    pmb->m0[i] = 0.;
    pmb->mevap[i] = 0.;
    pmb->mf[i] = 0.;
  }
  pmb->first = true;
  return pmb;
}

/**
### *delete_phase_mass_balances()*: free mass balance object */

void delete_phase_mass_balances (PhaseMassBalance * pmb) {
  free (pmb->m), pmb->m = NULL;
  free (pmb->m0), pmb->m = NULL;
  free (pmb->mevap), pmb->mevap = NULL;
  free (pmb->mf), pmb->mf = NULL;
  free (pmb), pmb = NULL;
}

/**
### *phase_mass_balance()*: compute mass balance for a single time step */

static double face_value_bid (Point point, scalar s, int bid) {
  double val = 0.;
  switch (bid) {
    case 0: val = 0.5*(s[1,0]  + s[]); break;   // right
    case 1: val = 0.5*(s[-1,0] + s[]); break;   // left
    case 2: val = 0.5*(s[0,1]  + s[]); break;   // top
    case 3: val = 0.5*(s[0,-1] + s[]); break;   // bottom
  }
  return val;
}

static double face_gradient_bid (Point point, scalar s, int bid) {
  double grad = 0.;
  switch (bid) {
    case 0: grad = (s[1,0]  - s[])/Delta*fm.x[];  break;  // right
    case 1: grad = (s[-1,0] - s[])/Delta*fm.x[];  break;  // left
    case 2: grad = (s[0,1]  - s[])/Delta*fm.y[];  break;  // top
    case 3: grad = (s[0,-1] - s[])/Delta*fm.y[];  break;  // bottom
  }
  return grad;
}

static double face_flux_bid (Point point, face vector uf, int bid) {
  double flux = 0.;
  switch (bid) {
    case 0: flux = +uf.x[]; break;   // right
    case 1: flux = -uf.x[]; break;   // left
    case 2: flux = +uf.y[]; break;   // top
    case 3: flux = -uf.y[]; break;   // bottom
  }
  return flux;
}

void phase_mass_balance (PhaseMassBalance * pmb, const Phase * phase,
    scalar * mEvapList = NULL, double dt, face vector uf,
    (const) scalar f = unity, bool boundaries = true,
    bool fick_corrected = false, bool molar_diffusion = false)
{
  PhaseMassBalance * balance = new_phase_mass_balance (phase);

  foreach_scalar_in (phase) {
    foreach(serial) {
      double ff = phase->inverse ? 1. - f[] : f[];
      balance->mtot += rho[]*ff*dv();
      foreach_species_in (phase)
        balance->m[i] += rho[]*Y[]*dv();
    }
    if (mEvapList) {  // fixme: must be computed befored the interface moves
      assert (phase->n == list_len (mEvapList));
      foreach_interfacial_plic (f, F_ERR, serial) {
        foreach_species_in (phase) {
          // note: both dirac and dv() multiply by cm
          scalar mEvap = mEvapList[i];
          balance->mevap[i] += mEvap[]*dirac*dt*dv()/cm[];
        }
      }
    }
  }

  if (boundaries) {
    foreach_scalar_in (phase) {
      // note: do not remove calls to boundaries
      // `foreach_boundary` does not trigger the automatic boundary conditions
      boundary ({rho,f,MW});
      boundary ((scalar *){uf});
      boundary (phase->YList);
      boundary (phase->XList);
      boundary (phase->DList);

      for (int b = 0; b < nboundary; b++) {
        foreach_boundary(b, serial) {
          // Convective flux
          double rhof = face_value_bid (point, rho, b);
          double uff = face_flux_bid (point, uf, b);
          double ff = face_value_bid (point, f, b);
          ff = phase->inverse ? 1. - ff : ff;
          balance->mftot += rhof*ff*uff*Delta*dt;
          foreach_species_in (phase)
            balance->mf[i] += rhof*ff*face_value_bid (point, Y, b)*uff*Delta*dt;

          // Diffusive flux
          // todo: need check metrics
          double jcorr = 0.;
          if (fick_corrected) {
            foreach_species_in (phase) {
              double Df = face_value_bid (point, D, b);
              double gradYn = 0.;
              if (molar_diffusion) {
                double MWf = face_value_bid (point, MW, b);
                gradYn = (MWf > 0.) ?
                  phase->MWs[i]/MWf*face_gradient_bid (point, X, b) : 0.;
              }
              else
                gradYn = face_gradient_bid (point, Y, b);
              jcorr += -rhof*Df*gradYn*Delta;
            }
          }

          foreach_species_in (phase) {
            double Df = face_value_bid (point, D, b);
            double Yf = face_value_bid (point, Y, b);
            double gradYn = 0.;
            if (molar_diffusion) {
              double MWf = face_value_bid (point, MW, b);
              gradYn = (MWf > 0.) ?
                phase->MWs[i]/MWf*face_gradient_bid (point, X, b) : 0.;
            }
            else
              gradYn = face_gradient_bid (point, Y, b);
            balance->mf[i] += (-rhof*Df*gradYn*Delta - jcorr*Yf)*dt;
          }
        }
      }
    }
  }

#if _MPI
  mpi_all_reduce (balance->mtot, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (balance->mftot, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce_array (balance->m, MPI_DOUBLE, MPI_SUM, phase->n);
  mpi_all_reduce_array (balance->mf, MPI_DOUBLE, MPI_SUM, phase->n);
  mpi_all_reduce_array (balance->mevap, MPI_DOUBLE, MPI_SUM, phase->n);
#endif

  if (pmb->first) {
    pmb->mtot0 = balance->mtot;
    foreach_species_in (phase)
      pmb->m0[i] = balance->m[i];
  }
  pmb->first = false;

  pmb->mtot = balance->mtot;
  pmb->mftot += balance->mftot;
  foreach_species_in (phase) {
    pmb->m[i] = balance->m[i];
    pmb->mevap[i] += balance->mevap[i];
    pmb->mf[i] += balance->mf[i];
    pmb->mevaptot += balance->mevap[i];
  }

  delete_phase_mass_balances (balance);
}

/**
## Notes

The C language does not allow inheritance. However, if we need to implement a
different type of phase with additional attributes, we can define an additional
structure:

~~~literatec
typedef struct {
  Phase phase;
  scalar new_attribute;
} CustomPhase;
~~~

Since the address of the structure is the same as the address of its first
element, the functions for the original `Phase` can be utilized also for the
`CustomPhase`, after proper recasting:

~~~literatec
CustomPhase * cphase = new_phase();
phase_function ((Phase *)cphase);
~~~
*/

