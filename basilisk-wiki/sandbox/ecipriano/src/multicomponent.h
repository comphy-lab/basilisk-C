#include "intgrad.h"
#include "fsolve.h"

Phase * liq_int, * gas_int;

scalar * mEvapList = NULL;

int * LSI = NULL, * GOSI = NULL, NGOS, inertIndex;

double * YIntVals;

typedef struct {
  Phase * liq, * gas;
  scalar fl, fg;
  face vector fsl, fsg;
  coord c;
} EnergyBalanceData;

double divq_rad_int (double TInti, double Tbulk = 300., double alphacorr = 0) {
  return alphacorr*5.670373e-8*(pow (Tbulk, 4.) - pow (TInti, 4.));
}

void energy_balance (const double * xdata, double * fdata, void * params) {
  EnergyBalanceData * data = params;

  Phase * liq = data->liq, * gas = data->gas;
  scalar fl = data->fl, fg = data->fg;
  face vector fsl = data->fsl, fsg = data->fsg;

  scalar TL = liq->T, TG = gas->T;

  foreach_point (data->c.x, data->c.y, data->c.z, serial) {
    double TInti = xdata[0];
    bool success = false;

    double gtrgrad = ebmgrad (point, TG, fl, fg, fsl, fsg, true,  TInti, &success);
    double ltrgrad = ebmgrad (point, TL, fl, fg, fsl, fsg, false, TInti, &success);

    scalar lambda1 = liq->lambda, lambda2 = gas->lambda;

    double vapheat = 0.;
    for (size_t i = 0; i < liq->n; i++) {
      scalar mEvap = mEvapList[LSI[i]];
      scalar dhev = liq->dhevList[i];
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

// Antoine functions and default condition

double antoine_default (double T, double P, int i) {
  if (YIntVal)
    return YIntVal;
  else {
    scalar YL = liq->YList[i];
    if (YL.antoine)
      return YL.antoine (T, P);
    else
      return YIntVals[i];
  }
}

double (* antoine) (double, double, int) = &antoine_default;

event defaults (i = 0) {
  liq_int = new_phase_minimal ("LInt", NLS, false, liq_species);
  gas_int = new_phase_minimal ("GInt", NGS, true, gas_species);

  foreach_species_in (gas_int) {
    scalar mEvap = new scalar;
    char name[80];
    sprintf (name, "mEvap%s", gas_int->species[i]);
    free (mEvap.name);
    mEvap.name = strdup (name);
    mEvap.nodump = true;
    mEvapList = list_add (mEvapList, mEvap);
  }

  // Liquid species indices *LSI* array
  Array * arrLSI = array_new();
  for (int i = 0; i < NGS; i++)
    for (int j = 0; j < NLS; j++)
      if (strcmp (gas_int->species[i], liq_int->species[j]) == 0)
        array_append (arrLSI, &i, sizeof (int));
  LSI = (int *)array_shrink (arrLSI);

  // Gas-only species indices *GOSI* array
  Array * arrGOSI = array_new();
  for (int i = 0; i < NGS; i++) {
    bool species_is_also_liquid = false;
    for (int j = 0; j < NLS; j++) {
      if (strcmp (gas_int->species[i], liq_int->species[j]) == 0)
        species_is_also_liquid = true;
    }
    if (!species_is_also_liquid)
      array_append (arrGOSI, &i, sizeof (int));
  }
  NGOS = arrGOSI->len / sizeof (int);
  GOSI = (int *)array_shrink (arrGOSI);
  assert (NGOS == (NGS - NLS));

  // fixme: needs the inert_species name
  inertIndex = NGS - 1;
  //for (int i = 0; i < NGS; i++)
  //  if (strcmp (gas_int->species[i], inert_species) == 0)
  //    inertIndex = i;

  // Default antoine functions: NLS and the species
  YIntVals = (double *)calloc (NLS, sizeof (double));

  // We dump only non-null vaporization rates
  foreach_species_in (liq_int) {
    scalar mEvap = mEvapList[LSI[i]];
    mEvap.nodump = false;
  }
}

event cleanup (t = end) {
  delete_phase (liq_int), liq_int = NULL;
  delete_phase (gas_int), gas_int = NULL;
  delete (mEvapList), free (mEvapList), mEvapList = NULL;
  free (LSI), LSI = NULL;
  free (GOSI), GOSI = NULL;
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

  // Calculate interfacial temperature
  scalar TLInt = liq_int->T, TGInt = gas_int->T;
  foreach_interfacial (f, F_ERR) {
    TLInt[] = TIntVal;
    TGInt[] = TIntVal;
    if (!pcm.isothermal_interface) {
      TLInt[] = avg_neighbor (point, liq->T, f);
      TGInt[] = TLInt[];
    }
  }

  // Calculate interfacial liquid mass fraction
  foreach_interfacial (f, F_ERR) {
    foreach_species_in (liq_int)
      Y[] = avg_neighbor (point, liq->YList[i], f);
    foreach_species_in (gas_int) {
      //Y[] = avg_neighbor (point, gas->YList[i], f);
      scalar YG = gas->YList[i];
      Y[] = YG[];
    }
  }

  // Update mole fractions and moles
  memcpy (liq_int->MWs, liq->MWs, liq->n*sizeof (double));
  memcpy (gas_int->MWs, gas->MWs, gas->n*sizeof (double));
  phase_update_mw_moles (liq_int);
  phase_update_mw_moles (gas_int);

  // Variables for interfacial mole fractions
  double * ygcorr = malloc ((NLS + 1)*sizeof (double));
  double * xgcorr = malloc ((NLS + 1)*sizeof (double));
  double * mwcorr = malloc ((NLS + 1)*sizeof (double));

  // Variables for gas-only mass fraction
  double * ygolist  = malloc (NGOS*sizeof (double));
  double * mwgolist = malloc (NGOS*sizeof (double));

  foreach_interfacial (f, F_ERR, serial) {

    // Store the sum of the gas-only species
    double sum_ngos_yg_old = 0., sum_ngos_xg_old = 0.;
    for (int i = 0; i < NGOS; i++) {
      scalar YGInt = gas_int->YList[GOSI[i]];
      scalar XGInt = gas_int->XList[GOSI[i]];
      sum_ngos_yg_old += YGInt[];
      sum_ngos_xg_old += XGInt[];
    }

    for (int i = 0; i < NLS; i++) {
      scalar XLInt = liq_int->XList[i];
      scalar XGInt = gas_int->XList[LSI[i]];

      XGInt[] = min (antoine (TLInt[], Pref, i), 0.98)*XLInt[];
    }

    // Calculate the molecular weight of the gas-only species mixture
    for (int i = 0; i < NGOS; i++) {
      scalar YGInt = gas_int->YList[GOSI[i]];
      ygolist[i] = YGInt[];
      mwgolist[i] = gas_int->MWs[GOSI[i]];
    }
    correctfrac (ygolist, NGOS);
    double mwgasonly = mass2mw (ygolist, mwgolist, NGOS);

    // Convert interfacial mole fractions to mass fractions
    double sum_xgl = 0.;
    for (int i = 0; i < NLS; i++) {
      scalar XGInt = gas_int->XList[LSI[i]];
      sum_xgl += XGInt[];
      xgcorr[i] = XGInt[];
      mwcorr[i] = liq_int->MWs[i];
    }
    xgcorr[NLS] = 1. - sum_xgl;
    mwcorr[NLS] = mwgasonly;
    mole2massfrac (ygcorr, xgcorr, mwcorr, NLS + 1);

    // Recover the interface gas-phase side mass fractions
    for (int i = 0; i < NLS; i++) {
      scalar YGInt = gas_int->YList[LSI[i]];
      YGInt[] = ygcorr[i];
    }

    // Adjust the interface gas-phase (gas only) mass fractions
    for (int i = 0; i < NGOS; i++) {
      scalar YGInt = gas_int->YList[GOSI[i]];
      scalar XGInt = gas_int->XList[GOSI[i]];
      YGInt[] *= ygcorr[NLS] / (sum_ngos_yg_old + 1e-10);
      XGInt[] *= xgcorr[NLS] / (sum_ngos_xg_old + 1e-10);
    }
  }
  free (ygcorr);
  free (xgcorr);
  free (mwcorr);
  free (ygolist);
  free (mwgolist);

  // Compute the total and species vaporization rates
  foreach_interfacial (f, F_ERR) {
    scalar rhoG = gas->rho;

    double jGtot = 0.;
    if (pcm.fick_corrected) {
      for (int i = 0; i < NGS; i++) {
        scalar DG = gas->DList[i];
        double gtrgrad = 0.;
        if (pcm.molar_diffusion) {
          scalar XG = gas->XList[i];
          scalar XGInt = gas_int->XList[i];
          scalar MWGInt = gas_int->MW;
          gtrgrad = ebmgrad (point, XG, fl, fg, fsl, fsg, true, XGInt[], false);
          gtrgrad *= (MWGInt[] > 0) ? gas_int->MWs[i]/MWGInt[] : 0.;
        }
        else {
          scalar YGInt = gas_int->YList[i];
          scalar YG = gas->YList[i];
          gtrgrad = ebmgrad (point, YG, fl, fg, fsl, fsg, true, YGInt[], false);
        }
        jGtot += -rhoG[]*DG[]*gtrgrad;
      }
    }

    double sum_jG = 0., sum_YGInt = 0.;
    for (int i = 0; i < NLS; i++) {
      scalar YGInt = gas_int->YList[LSI[i]];
      scalar DG = gas->DList[LSI[i]];
      double gtrgrad = 0.;
      if (pcm.molar_diffusion) {
        scalar XGInt = gas_int->XList[LSI[i]];
        scalar XG = gas->XList[LSI[i]];
        scalar MWGInt = gas_int->MW;

        gtrgrad = ebmgrad (point, XG, fl, fg, fsl, fsg, true, XGInt[], false);
        gtrgrad *= (MWGInt[] > 0) ? gas_int->MWs[LSI[i]]/MWGInt[] : 0.;
      }
      else {
        scalar YG = gas->YList[LSI[i]];
        gtrgrad = ebmgrad (point, YG, fl, fg, fsl, fsg, true, YGInt[], false);
      }
      sum_jG += -rhoG[]*DG[]*gtrgrad - YGInt[]*jGtot;
      sum_YGInt += YGInt[];
    }
    mEvapTot[] = sum_jG / min (1. - sum_YGInt, 0.99);

    for (int i = 0; i < NLS; i++) {
      scalar mEvap = mEvapList[LSI[i]];
      scalar YGInt = gas_int->YList[LSI[i]];
      scalar DG = gas->DList[LSI[i]];
      double gtrgrad = 0.;
      if (pcm.molar_diffusion) {
        scalar XGInt = gas_int->XList[LSI[i]];
        scalar XG = gas->XList[LSI[i]];
        scalar MWGInt = gas_int->MW;

        gtrgrad = ebmgrad (point, XG, fl, fg, fsl, fsg, true, XGInt[], false);
        gtrgrad *= (MWGInt[] > 0) ? gas_int->MWs[LSI[i]]/MWGInt[] : 0.;
      }
      else {
        scalar YG = gas->YList[LSI[i]];
        gtrgrad = ebmgrad (point, YG, fl, fg, fsl, fsg, true, YGInt[], false);
      }
      mEvap[] = mEvapTot[]*YGInt[] - rhoG[]*DG[]*gtrgrad - YGInt[]*jGtot;
    }
  }

  // Find the interfacial temperature which respects the vaporization rate
  if (!pcm.isothermal && !pcm.isothermal_interface) {
    EnergyBalanceData data;
    data.liq = liq, data.gas = gas;
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
          "src/multicomponent.h:%d: error: missing root-finding method\n",
          LINENO), fflush (ferr);
      exit(1);
#endif

      double * unk = (double *)arrUnk->p;
      TLInt[] = unk[0];
      TGInt[] = TLInt[];
      array_free (arrUnk);
    }
  }

  // Calculate species and temperature source terms
  foreach_interfacial_plic (f, F_ERR) {

    foreach_species_in (gas) {
      scalar mEvap = mEvapList[i];
      scalar YInt = gas_int->YList[i];

      if (pcm.diffusion == EXPLICIT_ONLY)
        SYexp[] += -(mEvap[] - mEvapTot[]*YInt[])*dirac;
      else {
        SYexp[] += -mEvap[]*dirac;
        SYimp[] += +mEvapTot[]*dirac;
      }
    }

    foreach_species_in (liq) {
      scalar mEvap = mEvapList[LSI[i]];
      scalar YInt = liq_int->YList[i];

      if (pcm.diffusion == EXPLICIT_ONLY)
        SYexp[] += +(mEvap[] - mEvapTot[]*YInt[])*dirac;
      else {
        SYexp[] += +mEvap[]*dirac;
        SYimp[] += -mEvapTot[]*dirac;
      }
    }

    scalar TL = liq->T, TG = gas->T;
    double ltrgrad = ebmgrad (point, TL, fl, fg, fsl, fsg, false, TLInt[], false);
    double gtrgrad = ebmgrad (point, TG, fl, fg, fsl, fsg, true, TGInt[], false);

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

