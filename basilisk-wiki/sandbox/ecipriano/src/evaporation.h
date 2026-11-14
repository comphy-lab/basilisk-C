/**
# Evaporation

Diffusion-driven phase change of multicomponent mixtures in non-isothermal
conditions, following [Cipriano et al., 2023](#cipriano2023), and [Cipriano
et al., 2024](#cipriano2024) for variable properties. This module
includes [phasechange.h](phasechange.h), computes the vaporization rates, and
sets the phase change model parameters (`pcm`) for evaporation. In each
interfacial cell, it solves the interface jump condition, assuming
thermodynamic equilibrium at the interface:
$$
  \begin{cases}
    \dot{m}_i =
    \dot{m}\hat\omega_{i,l}
    + \rho_l\mathcal{D}_{i,l}\nabla\omega_{i,l}\cdot\mathbf{n}_\Gamma
    =
    \dot{m}\hat\omega_{i,g}
    + \rho_g\mathcal{D}_{i,g}\nabla\omega_{i,g}\cdot\mathbf{n}_\Gamma
    \\
    \sum_i \dot{m}_i\Delta h_{ev,i} =
    \lambda_l\nabla T_l\cdot\mathbf{n}_\Gamma
    - \lambda_g\nabla T_g\cdot\mathbf{n}_\Gamma
    \\
    \hat{x}_{i,g} = K_{eq,i}(\hat{T})\,\hat{x}_{i,l}
  \end{cases}
$$
where the hat denotes interface values, and $\mathbf{n}_\Gamma$ is the
interface normal pointing towards the liquid. The equilibrium is written here
for evaporation; in general, it relates the continuous to the dispersed mole
fractions (see below). The solution gives the vaporization rate
$\dot{m}_i$ of each species, the interface mass fractions
$\hat\omega_{i,l}$ and $\hat\omega_{i,g}$, and the interface temperature
$\hat{T}$, which provide the source terms of the species and temperature
equations solved by [phasechange.h](phasechange.h).

## Dispersed and continuous phases

The two phases play different roles in the phase change:

* the *dispersed* phase contains only species that cross the interface;
* the *continuous* phase contains the same species, plus *inert* species that
  never cross the interface.

In droplet evaporation, the liquid fuel droplet is the dispersed phase, and the
gas (fuel vapor plus air, the inert) is the continuous phase. In bubble
dissolution, or bubble growth during electrolysis, the gas bubble is the
dispersed phase, and the liquid (dissolved gas plus water, the inert) is the
continuous phase. The names do not imply a direction: the mass flux is
computed from the dispersed to the continuous phase, and it is negative if
species move the other way, as in a growing bubble.

The dispersed phase is identified as the phase with fewer species (the liquid
if `NLS` and `NGS` are equal), not from the topology of the interface: in
planar or stratified configurations (e.g. the Stefan problem), it is simply
the phase without inert species. Each dispersed species must also be a
continuous species. The following variables describe the two sets of species:

* *NDS*: number of dispersed species (all the species of the dispersed phase);
* *NCS*: number of continuous species;
* *NIS*: number of inert species, *NCS - NDS*;
* *DSI*: dispersed species indices, i.e. the continuous index of each
  dispersed species;
* *ISI*: inert species indices, i.e. the continuous index of each inert
  species.

The roles are chosen once, in *liquid_dispersed*. The species part of the jump
condition is computed by *species_jump_condition()*, which takes the dispersed
and the continuous phases as arguments: it is called with the liquid as the
dispersed phase for evaporation, and with the gas as the dispersed phase for
bubble dissolution.

The vaporization rates `mEvapList` are indexed like the continuous species. The
sign convention is the one in [phasechange.h](phasechange.h): `mEvapTot` and
`mEvapList` are positive when the liquid loses mass. */

#include "phasechange.h"
#include "gradients.h"
#include "fsolve.h"

Phase * liq_int, * gas_int;

scalar * mEvapList = NULL;

bool liquid_dispersed = true;
int NDS = 0, NCS = 0, NIS = 0;
int * DSI = NULL, * ISI = NULL;

double * YIntVals;

/**
## Interface energy balance */

typedef struct {
  Phase * liq, * gas, * dispersed;
  scalar fl, fg;
  face vector fsl, fsg;
  coord c;
} EnergyBalanceData;

double divq_rad_int (double TInti, double Tbulk = 300., double alphacorr = 0) {
  return alphacorr*5.670373e-8*(pow (Tbulk, 4.) - pow (TInti, 4.));
}

void energy_balance (const double * xdata, double * fdata, void * params) {
  EnergyBalanceData * data = params;

  Phase * liq = data->liq, * gas = data->gas, * d = data->dispersed;
  scalar fl = data->fl, fg = data->fg;
  face vector fsl = data->fsl, fsg = data->fsg;

  scalar TL = liq->T, TG = gas->T;

  foreach_point (data->c.x, data->c.y, data->c.z, serial) {
    double TInti = xdata[0];

    double gtrgrad = plic_gradient (point, TG, fg, fsg, TInti, true, NULL);
    double ltrgrad = plic_gradient (point, TL, fl, fsl, TInti, true, NULL);

    scalar lambda1 = liq->lambda, lambda2 = gas->lambda;

    double vapheat = 0.;
    for (size_t i = 0; i < d->n; i++) {
      scalar mEvap = mEvapList[DSI[i]];
      scalar dhev = d->dhevList[i];
      vapheat -= mEvap[]*dhev[];
    }

    fdata[0] = vapheat
      - divq_rad_int (TInti, TG0, pcm.emissivity)
      + lambda1[]*ltrgrad
      + lambda2[]*gtrgrad
      ;
  }
}

#if USE_GSL
int energy_balance_gsl (const gsl_vector * x, void * params, gsl_vector * f) {
  double * xdata = x->data, * fdata = f->data;
  energy_balance (xdata, fdata, params);
  return GSL_SUCCESS;
}
#endif

/**
## Interface equilibrium

The function *antoine()* returns the equilibrium constant of the i-th
dispersed species, i.e. the ratio between its continuous-side and
dispersed-side interface mole fractions. For evaporation it is Raoult's law,
$p^{sat}_i(T)/P$; for bubble dissolution it is Henry's law, $P/H_i(T)$. The
default implementation uses the *antoine* attribute of the dispersed mass
fraction, if set. */

double antoine_default (double T, double P, int i) {
  if (YIntVal)
    return YIntVal;
  else {
    Phase * d = liquid_dispersed ? liq : gas;
    scalar YD = d->YList[i];
    if (YD.antoine)
      return YD.antoine (T, P);
    else
      return YIntVals[i];
  }
}

double (* antoine) (double, double, int) = &antoine_default;

/**
## Species jump condition

Normal gradient, at the interface, of the continuous species *k*, from mass or
mole fractions. */

double continuous_gradient (Point point, Phase * c, Phase * c_int, int k,
    scalar fc, face vector fsc)
{
  double grad = 0.;
  if (pcm.molar_diffusion) {
    scalar X = c->XList[k], XInt = c_int->XList[k], MWInt = c_int->MW;
    grad = plic_gradient (point, X, fc, fsc, XInt[], true, NULL);
    grad *= (MWInt[] > 0) ? c_int->MWs[k]/MWInt[] : 0.;
  }
  else {
    scalar Y = c->YList[k], YInt = c_int->YList[k];
    grad = plic_gradient (point, Y, fc, fsc, YInt[], true, NULL);
  }
  return grad;
}

/**

The function below solves the species part of the jump condition, from the
dispersed phase *d* to the continuous phase *c*, with the decoupled procedure
of [Cipriano et al., 2023](#cipriano2023), Sec. 3.5: the dispersed phase is
assumed to evolve more slowly than the continuous phase (e.g. the liquid
during evaporation), so its interface values are taken from the surrounding
cells, while the continuous side follows from equilibrium and mass balances.
*d_int* and *c_int* store the interface values, *fc* and *fsc* are the volume
and face fractions of the continuous phase, and *sign* converts the
dispersed-to-continuous mass flux into vaporization rates (+1 if the dispersed
phase is the liquid). The interface temperature must be already set in
*d_int*. */

void species_jump_condition (Phase * d, Phase * c, Phase * d_int,
    Phase * c_int, scalar fc, face vector fsc, double sign)
{

  /**
  The dispersed interface mass fractions are the volume-fraction-weighted
  averages over the 3x3 stencil (Eq. 45). The continuous values only provide
  the inert composition, rescaled below. */

  foreach_interfacial (f, F_ERR) {
    foreach_species_in (d_int)
      Y[] = avg_neighbor (point, d->YList[i], f);
    foreach_species_in (c_int) {
      scalar YC = c->YList[i];
      Y[] = YC[];
    }
  }

  // Update mole fractions and moles
  memcpy (d_int->MWs, d->MWs, d->n*sizeof (double));
  memcpy (c_int->MWs, c->MWs, c->n*sizeof (double));
  phase_update_mw_moles (d_int);
  phase_update_mw_moles (c_int);

  /**
  The continuous interface mole fractions follow from the phase equilibrium
  (Eq. 47)
  $$
  \hat{x}_{i,c} = K_{eq,i}(\hat{T})\,\hat{x}_{i,d}
  $$
  with $K_{eq,i}$ given by *antoine()*. They are converted to mass fractions
  by lumping the inert species into a single pseudo-species (the last element
  of the arrays below), which closes the mole fractions to one. The inert
  species are then rescaled, preserving their relative composition. */

  scalar TInt = d_int->T;
  foreach_interfacial (f, F_ERR, serial) {
    double yccorr[NDS + 1], xccorr[NDS + 1], mwcorr[NDS + 1];
    double yinert[NIS > 0 ? NIS : 1], mwinert[NIS > 0 ? NIS : 1];

    // Store the sum of the inert species
    double sum_yinert_old = 0., sum_xinert_old = 0.;
    for (int j = 0; j < NIS; j++) {
      scalar YCInt = c_int->YList[ISI[j]];
      scalar XCInt = c_int->XList[ISI[j]];
      sum_yinert_old += YCInt[];
      sum_xinert_old += XCInt[];
    }

    // Thermodynamic equilibrium
    for (int i = 0; i < NDS; i++) {
      scalar XDInt = d_int->XList[i];
      scalar XCInt = c_int->XList[DSI[i]];
      XCInt[] = min (antoine (TInt[], P0, i), 0.98)*XDInt[];
    }

    // Molecular weight of the inert species mixture
    for (int j = 0; j < NIS; j++) {
      scalar YCInt = c_int->YList[ISI[j]];
      yinert[j] = YCInt[];
      mwinert[j] = c_int->MWs[ISI[j]];
    }
    correctfrac (yinert, NIS);
    double mwinertmix = mass2mw (yinert, mwinert, NIS);

    // Convert interface mole fractions to mass fractions
    double sum_xc = 0.;
    for (int i = 0; i < NDS; i++) {
      scalar XCInt = c_int->XList[DSI[i]];
      sum_xc += XCInt[];
      xccorr[i] = XCInt[];
      mwcorr[i] = d_int->MWs[i];
    }
    xccorr[NDS] = 1. - sum_xc;
    mwcorr[NDS] = mwinertmix;
    mole2massfrac (yccorr, xccorr, mwcorr, NDS + 1);

    for (int i = 0; i < NDS; i++) {
      scalar YCInt = c_int->YList[DSI[i]];
      YCInt[] = yccorr[i];
    }

    // Rescale the inert species
    for (int j = 0; j < NIS; j++) {
      scalar YCInt = c_int->YList[ISI[j]];
      scalar XCInt = c_int->XList[ISI[j]];
      YCInt[] *= yccorr[NDS] / (sum_yinert_old + 1e-10);
      XCInt[] *= xccorr[NDS] / (sum_xinert_old + 1e-10);
    }
  }

  /**
  Since the dispersed diffusive fluxes sum to zero, summing the species
  balances over the dispersed species gives an explicit expression for the
  total mass flux from the dispersed to the continuous phase (Eqs. 49, 50)
  $$
  \dot{m} = \frac{\sum_{i=1}^{NDS} j_{i,c}}
    {1 - \sum_{i=1}^{NDS} \hat{\omega}_{i,c}}
  \quad
  \dot{m}_i = \dot{m}\,\hat{\omega}_{i,c} + j_{i,c}
  $$
  where $j_{i,c}$ is the diffusive flux of the continuous species along the
  normal $\mathbf{n}_c$ entering the continuous phase. It is defined by two
  options ([Cipriano et al., 2024](#cipriano2024), Eqs. 26, 27):

  * `pcm.molar_diffusion`: the diffusivities are molar based (e.g. from the
  Hirschfelder and Curtiss approximation), so Fick's law is written in terms of
  mole fractions, $j^F_{i,c} = -\rho_c\mathcal{D}_{i,c} (M_{w,i}/M_{w,c})\,
  \partial x_{i,c}/\partial\mathbf{n}_c$, with $M_{w,c}$ the molecular weight
  of the interface mixture. Otherwise, $j^F_{i,c} = -\rho_c\mathcal{D}_{i,c}\,
  \partial\omega_{i,c}/\partial\mathbf{n}_c$.
  * `pcm.fick_corrected`: Fick's law is exact only for binary mixtures, or for
  equal diffusivities. The correction
  $j_{i,c} = j^F_{i,c} - \hat{\omega}_{i,c}\sum_{k=1}^{NCS} j^F_{k,c}$
  makes the diffusive fluxes sum to zero. Otherwise, $j_{i,c} = j^F_{i,c}$.
  */

  foreach_interfacial (f, F_ERR) {
    scalar rhoC = c->rho;

    double jCtot = 0.;
    if (pcm.fick_corrected) {
      for (int k = 0; k < NCS; k++) {
        scalar DC = c->DList[k];
        double ctrgrad = continuous_gradient (point, c, c_int, k, fc, fsc);
        jCtot += -rhoC[]*DC[]*ctrgrad;
      }
    }

    double sum_jC = 0., sum_YCInt = 0.;
    for (int i = 0; i < NDS; i++) {
      scalar YCInt = c_int->YList[DSI[i]];
      scalar DC = c->DList[DSI[i]];
      double ctrgrad = continuous_gradient (point, c, c_int, DSI[i], fc, fsc);
      sum_jC += -rhoC[]*DC[]*ctrgrad - YCInt[]*jCtot;
      sum_YCInt += YCInt[];
    }
    double mtot = sum_jC / min (1. - sum_YCInt, 0.99);
    mEvapTot[] = sign*mtot;

    for (int i = 0; i < NDS; i++) {
      scalar mEvap = mEvapList[DSI[i]];
      scalar YCInt = c_int->YList[DSI[i]];
      scalar DC = c->DList[DSI[i]];
      double ctrgrad = continuous_gradient (point, c, c_int, DSI[i], fc, fsc);
      mEvap[] = sign*(mtot*YCInt[] - rhoC[]*DC[]*ctrgrad - YCInt[]*jCtot);
    }
  }

  /**
  The species equations get the explicit source $\dot{m}_i\delta_\Gamma$ and
  the implicit source $\dot{m}\,\omega_i\delta_\Gamma$, which accounts for
  the volume change (Eq. 43; both explicit with `EXPLICIT_ONLY`). Since
  vaporization rates are positive when the liquid loses mass, the sign of the
  sources is *sign* for the dispersed phase and *-sign* for the continuous
  phase. */

  foreach_interfacial_plic (f, F_ERR) {
    for (int k = 0; k < NCS; k++) {
      scalar mEvap = mEvapList[k], YInt = c_int->YList[k];
      scalar SYexp = c->SYexpList[k], SYimp = c->SYimpList[k];
      if (pcm.diffusion == EXPLICIT_ONLY)
        SYexp[] += -sign*(mEvap[] - mEvapTot[]*YInt[])*dirac;
      else {
        SYexp[] += -sign*mEvap[]*dirac;
        SYimp[] += sign*mEvapTot[]*dirac;
      }
    }

    for (int i = 0; i < NDS; i++) {
      scalar mEvap = mEvapList[DSI[i]], YInt = d_int->YList[i];
      scalar SYexp = d->SYexpList[i], SYimp = d->SYimpList[i];
      if (pcm.diffusion == EXPLICIT_ONLY)
        SYexp[] += sign*(mEvap[] - mEvapTot[]*YInt[])*dirac;
      else {
        SYexp[] += sign*mEvap[]*dirac;
        SYimp[] += -sign*mEvapTot[]*dirac;
      }
    }
  }
}

/**
## Setup */

event defaults (i = 0) {
  liq_int = new_phase_minimal ("LInt", NLS, false, liq_species);
  gas_int = new_phase_minimal ("GInt", NGS, true, gas_species);

  liquid_dispersed = (NLS <= NGS);
  Phase * d_int = liquid_dispersed ? liq_int : gas_int;
  Phase * c_int = liquid_dispersed ? gas_int : liq_int;

  NDS = d_int->n;
  NCS = c_int->n;
  NIS = NCS - NDS;

  foreach_species_in (c_int) {
    scalar mEvap = new scalar;
    char name[80];
    sprintf (name, "mEvap%s", c_int->species[i]);
    free (mEvap.name);
    mEvap.name = strdup (name);
    mEvap.nodump = true;
    mEvapList = list_add (mEvapList, mEvap);
  }

  // Dispersed species indices *DSI* (aborts if missing from the continuous
  // phase)
  DSI = (int *)malloc (NDS*sizeof (int));
  for (int i = 0; i < NDS; i++)
    DSI[i] = phase_species_index (c_int, d_int->species[i]);

  // Inert species indices *ISI*: the continuous species not in *DSI*
  bool dispersed_species[NCS];
  for (int k = 0; k < NCS; k++)
    dispersed_species[k] = false;
  for (int i = 0; i < NDS; i++)
    dispersed_species[DSI[i]] = true;

  ISI = (int *)malloc (NIS*sizeof (int));
  int nis = 0;
  for (int k = 0; k < NCS; k++)
    if (!dispersed_species[k])
      ISI[nis++] = k;
  assert (nis == NIS);

  // Default equilibrium constants of the dispersed species
  YIntVals = (double *)calloc (NDS, sizeof (double));

  // We dump only non-null vaporization rates
  for (int i = 0; i < NDS; i++) {
    scalar mEvap = mEvapList[DSI[i]];
    mEvap.nodump = false;
  }
}

event cleanup (t = end) {
  delete_phase (liq_int), liq_int = NULL;
  delete_phase (gas_int), gas_int = NULL;
  delete (mEvapList), free (mEvapList), mEvapList = NULL;
  free (DSI), DSI = NULL;
  free (ISI), ISI = NULL;
  free (YIntVals);
}

event reset_sources (i++) {
  foreach() {
    foreach_scalar_in (liq_int) {
      T[] = 0.;
      foreach_species_in (liq_int) {
        Y[] = 0.;
        X[] = 0.;
      }
    }
    foreach_scalar_in (gas_int) {
      T[] = 0.;
      foreach_species_in (gas_int) {
        Y[] = 0.;
        X[] = 0.;
      }
    }
  }
}

/**
## Phase change */

event phasechange (i++) {
  phase_tracers_to_scalars (liq, f, tol = F_ERR);
  phase_tracers_to_scalars (gas, f, tol = F_ERR);

  scalar fl[], fg[];
  foreach() {
    fl[] = f[];
    fg[] = 1. - f[];
  }

  face vector fsl[], fsg[];
  face_fraction (fl, fsl);
  face_fraction (fg, fsg);

  /**
  The interface temperature is first estimated from the liquid temperature
  (or set to *TIntVal* if `pcm.isothermal_interface`), using the
  volume-fraction-weighted average over the 3x3 stencil. */

  scalar TLInt = liq_int->T, TGInt = gas_int->T;
  foreach_interfacial (f, F_ERR) {
    TLInt[] = TIntVal;
    TGInt[] = TIntVal;
    if (!pcm.isothermal_interface) {
      TLInt[] = avg_neighbor (point, liq->T, f);
      TGInt[] = TLInt[];
    }
  }

  if (liquid_dispersed)
    species_jump_condition (liq, gas, liq_int, gas_int, fg, fsg, 1.);
  else
    species_jump_condition (gas, liq, gas_int, liq_int, fl, fsl, -1.);

  /**
  Unless the interface is isothermal, the interface temperature is refined
  by solving the energy balance ([Cipriano et al., 2023](#cipriano2023),
  Eq. 51)
  $$
  \sum_{i=1}^{NDS} \dot{m}_i\Delta h_{ev,i} =
    \lambda_l\nabla T_l\cdot\mathbf{n}_\Gamma
    - \lambda_g\nabla T_g\cdot\mathbf{n}_\Gamma
  $$
  with $\mathbf{n}_\Gamma$ pointing towards the liquid, the vaporization
  rates frozen, and $\Delta h_{ev,i}$ taken from the dispersed species. In
  *energy_balance()*, each gradient is computed along the normal entering its
  own phase, hence the plus sign. */

  if (!pcm.isothermal && !pcm.isothermal_interface) {
    EnergyBalanceData data;
    data.liq = liq, data.gas = gas;
    data.dispersed = liquid_dispersed ? liq : gas;
    data.fl = fl, data.fg = fg;
    data.fsl = fsl, data.fsg = fsg;

    foreach_interfacial (f, F_ERR, serial) {
      Array * arrUnk = array_new();
      double vali = TLInt[];
      array_append (arrUnk, &vali, sizeof (double));

      coord o = {x, y, z};
      foreach_dimension()
        data.c.x = o.x;

#if USE_GSL
      fsolve (energy_balance_gsl, arrUnk, &data);
#else
      fprintf (ferr,
          "src/evaporation.h:%d: error: missing root-finding method\n",
          LINENO), fflush (ferr);
      exit(1);
#endif

      double * unk = (double *)arrUnk->p;
      TLInt[] = unk[0];
      TGInt[] = TLInt[];
      array_free (arrUnk);
    }
  }

  /**
  The temperature equations get the interface conductive fluxes ([Cipriano
  et al., 2023](#cipriano2023), Eq. 44). */

  foreach_interfacial_plic (f, F_ERR) {
    scalar TL = liq->T, TG = gas->T;
    double ltrgrad = plic_gradient (point, TL, fl, fsl, TLInt[], true, NULL);
    double gtrgrad = plic_gradient (point, TG, fg, fsg, TGInt[], true, NULL);

    scalar slT = liq->STexp, sgT = gas->STexp;
    scalar lambdal = liq->lambda, lambdag = gas->lambda;

    slT[] += lambdal[]*ltrgrad*dirac;
    sgT[] += lambdag[]*gtrgrad*dirac;
  }

  // We add the heat from the mass diffusion process
  if (pcm.mass_diffusion_enthalpy) {
    phase_add_heat_species_diffusion (liq, f, pcm.molar_diffusion, F_ERR);
    phase_add_heat_species_diffusion (gas, f, pcm.molar_diffusion, F_ERR);
  }

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);
}

/**
## Phase change model parameters

We set the `pcm` parameters for evaporation: the expansion term is shifted to
the liquid, and the species source terms are split into explicit and implicit
contributions. With a single velocity field (`nv == 1`), the interface
regression velocity uses the gas density. */

// fixme: check compatibility with restore
event init (i = 0) {
  pcm.shifting = SHIFT_TO_LIQUID;
  pcm.diffusion = EXPLICIT_IMPLICIT;
  if (nv == 1)
    pcm.byrhogas = true;
}

/**
## References

~~~bib
@hal{cipriano2023, hal-04252476}

@hal{cipriano2024, hal-04580590}
~~~
*/
