
#include "intgrad.h"
#include "fsolve-gsl.h"
#include "thermodynamics.h"

extern scalar * mEvapList;
extern scalar * YLList;
extern scalar * YGList;
extern scalar * YLIntList;
extern scalar * YGIntList;
#ifdef VARPROP
extern scalar * Dmix1List;
extern scalar * Dmix2List;
extern scalar * dhevList;
#endif
extern double inDmix1[NLS];
extern double inDmix2[NGS];
extern double inKeq[NLS];
extern double inMW[NGS];
extern double Pref;
extern scalar fL, fG;
extern face vector fsL, fsG;
extern int * LSI;
extern int * GOSI;
extern int inertIndex;
extern int NGOS;
extern double rho1, rho2;
#ifdef SOLVE_TEMPERATURE
extern double lambda1, lambda2, dhev, cp1, cp2, TL0, TG0;
extern scalar TInt, TL, TG, T;
#endif

#ifdef CLAPEYRON
extern double Tsat;
#endif
#ifdef USE_CLAPEYRON
extern double Tboil[NLS];
#endif

typedef struct {
  coord c;
} UserDataNls;

void Equations (const double * xdata, double * fdata, void * params) {
  UserDataNls * data = (UserDataNls *)params;

  double mEvapToti[NGS];
  double mEvapi[NLS];
  double YLInti[NLS];
  double YGInti[NGS];
  double XLInti[NLS];
  double XGInti[NGS];
  double jL[NLS];
  double jG[NGS];
  double mEvapSum;
  bool success;
  int count;
#ifdef SOLVE_TEMPERATURE
  double TInti;
  double gradTGn;
  double gradTLn;
#endif

  foreach_point (data->c.x, data->c.y, data->c.z, serial) {

    /**
    Rename unknowns for convenience. */

    count = 0;

    for (int jj=0; jj<NLS; jj++)
      mEvapi[jj] = xdata[count++];

    for (int jj=0; jj<NLS; jj++)
      YLInti[jj] = xdata[count++];

    for (int jj=0; jj<NGS; jj++)
      YGInti[jj] = xdata[count++];

#ifdef SOLVE_TEMPERATURE
    TInti = xdata[count++];
#endif

    /**
    Convert mass-to-mole fractions. */

    {
      double inMWG[NGS];
      for (int jj=0; jj<NGS; jj++)
        inMWG[jj] = inMW[jj];

      double inMWL[NLS];
      for (int jj=0; jj<NLS; jj++)
        inMWL[jj] = inMW[LSI[jj]];

      mass2molefrac (XLInti, YLInti, inMWL, NLS);
      mass2molefrac (XGInti, YGInti, inMWG, NGS);
    }

    /**
    Set to zero gas-only species vaporization mass flow-rates. */

    for (int jj=0; jj<NLS; jj++)
      mEvapToti[LSI[jj]] = mEvapi[jj];

    for (int jj=0; jj<NGOS; jj++)
      mEvapToti[GOSI[jj]] = 0.;

    /**
    Compute diffusive fluxes. */

    for (int jj=0; jj<NGS; jj++) {
      double rho2vh = rho2;
      double Dmix2vh = inDmix2[jj];
#ifdef VARPROP
      rho2vh = rho2v[];
      scalar Dmix2 = Dmix2List[jj];
      Dmix2vh = Dmix2[];
      //rho2vh = tp2.rhov (&ts2);
      //Dmix2vh = tp2.diff (&ts2, jj);
#endif
      scalar YG = YGList[jj];
      double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YGInti[jj], &success);
      jG[jj] = -rho2vh*Dmix2vh*gtrgrad;
    }
    for (int jj=0; jj<NLS; jj++) {
      double rho1vh = rho1;
      double Dmix1vh = inDmix1[jj];
#ifdef VARPROP
      rho1vh = rho1v[];
      scalar Dmix1 = Dmix1List[jj];
      Dmix1vh = Dmix1[];
      //rho1vh = tp1.rhov (&ts1);
      //Dmix1vh = tp1.diff (&ts1, jj);
#endif
      scalar YL = YLList[jj];
      double ltrgrad = ebmgrad (point, YL, fL, fG, fsL, fsG, false, YLInti[jj], &success);
      jL[jj] = rho1vh*Dmix1vh*ltrgrad;
    }

    /**
    Compute total diffusive fluxes for Fick corrected approach. */

    double jLtot = 0., jGtot = 0.;
#ifdef FICK_CORRECTED
    for (int jj=0; jj<NLS; jj++)
      jLtot += jL[jj];

    for (int jj=0; jj<NGS; jj++)
      jGtot += jG[jj];
#endif

#ifdef SOLVE_TEMPERATURE
    gradTGn = ebmgrad (point, TG, fL, fG, fsL, fsG, true,  TInti, &success);
    gradTLn = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TInti, &success);
#endif

    /**
    Compute sum of vaporization mass flow-rate. */

    mEvapSum = 0.;
    for (int jj=0; jj<NLS; jj++) {
      mEvapSum += mEvapi[jj];
    }

    count = 0;

    /**
    [NEW]
    Mass balance in gas phase for evaporating species. */

    for (int jj=0; jj<NLS; jj++) {

      fdata[count++] = mEvapi[jj]
                     - mEvapSum*YGInti[LSI[jj]]
                     - jG[jj]
                     + jGtot*YGInti[LSI[jj]]
                     ;
    }

    /**
    [NEW]
    Mass balance in liquid phase for liquid species. */

    if (NLS > 1) {
      for (int jj=0; jj<NLS; jj++) {

        fdata[count++] = mEvapi[jj]
                       - mEvapSum*YLInti[jj]
                       - jL[jj]
                       + jLtot*YLInti[jj];
                       ;
      }
    }
    else {
      fdata[count++] = YLInti[0] - 1.;
    }


    ///**
    //Mass balance for species in liquid phase. */

    //if (NLS > 1) {
    //  for (int jj=0; jj<NLS; jj++) {

    //    fdata[count++] = mEvapi[jj]
    //                   - mEvapSum*YLInti[jj]
    //                   - jL[jj]
    //                   + jLtot*YLInti[jj]
    //                   ;
    //  }
    //}
    //else {
    //  fdata[count++] = YLInti[0] - 1.;
    //}

    ///**
    //Mass balance for species in gas phase and for gas-only species. */

    //for (int jj=0; jj<NGS; jj++) {

    //  fdata[count++] = mEvapToti[jj]
    //                 - mEvapSum*YGInti[jj]
    //                 - jG[jj]
    //                 + jGtot*YGInti[jj]
    //                 ;
    //}

    /**
    Thermodynamic (VLE) equilibrium at the interface. */

    double Keq[NLS];
    for (int jj=0; jj<NLS; jj++)
      Keq[jj] = inKeq[jj];

#ifdef USE_CLAPEYRON
    for (int jj=0; jj<NLS; jj++) {
      Keq[jj] = min (clapeyron ( min (TInti, Tboil[jj]-1), Tboil[jj], dhev, inMW[LSI[jj]]), 0.98);
    }
#endif
#ifdef USE_ANTOINE
    for (int jj=0; jj<NLS; jj++) {
      scalar YL = YLList[jj];
      Keq[jj] = min (YL.antoine (TInti, Pref), 0.98);
    }
#endif
#ifdef USE_ANTOINE_OPENSMOKE
    for (int jj=0; jj<NLS; jj++) {
      Keq[jj] = min (opensmoke_antoine (TInti, Pref, jj), 0.98);
    }
#endif

    for (int jj=0; jj<NLS; jj++) {
      fdata[count++] = XLInti[jj]*Keq[jj] - XGInti[LSI[jj]];
    }

    /**
    [NEW]
    Mass balance in gas phase for gas-only species. */

    for (int jj=0; jj<NGOS; jj++) {
      fdata[count++] = mEvapToti[GOSI[jj]]
                     - mEvapSum*YGInti[GOSI[jj]]
                     - jG[GOSI[jj]]
                     + jGtot*YGInti[GOSI[jj]]
                     ;
    }

#ifdef SOLVE_TEMPERATURE

    /**
    Interface energy balance. */

    double vapheat = 0.;
    for (int jj=0; jj<NLS; jj++) {
      double dhevvh = dhev;
#ifdef VARPROP
      scalar dhevjj = dhevList[jj];
      dhevvh = dhevjj[];
#endif
      vapheat -= dhevvh*mEvapi[jj];
    }

    fdata[count++] = vapheat
# ifdef VARPROP
                   + lambda1v[]*gradTLn
                   + lambda2v[]*gradTGn
# else
                   + lambda1*gradTLn
                   + lambda2*gradTGn
# endif
                   ;
#endif
  }
}

int EquationsGsl (const gsl_vector * x, void * params, gsl_vector * f) {
  double * xdata = x->data;
  double * fdata = f->data;

  Equations (xdata, fdata, params);
  return GSL_SUCCESS;
}

/**
A partire da *YLIntList*, *YGIntList*, *mEvapList* e delle proprietà
fisiche come densità e diffusività, bisogna calcolare la jump condition. */

void ijc_CoupledNls ()
{
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      /**
      Reset variables. */

      Array * arrUnk = array_new();

      // mInt field
      for (int jj=0; jj<NLS; jj++) {
        scalar mEvapi = mEvapList[LSI[jj]];
        double vali = mEvapi[];
        array_append (arrUnk, &vali, sizeof(double));
      }

      // YLInt field
      for (int jj=0; jj<NLS; jj++) {
        //scalar YLi = YLList[jj];
        //double vali = YLi[]/max(f[], I_TOL);
        scalar YLInti = YLIntList[jj];
        double vali = YLInti[];
        array_append (arrUnk, &vali, sizeof(double));
      }

      // YGInt field
      for (int jj=0; jj<NGS; jj++) {
        //scalar YGi = YGList[jj];
        //double vali = YGi[]/max(1. - f[], I_TOL);
        scalar YGInti = YGIntList[jj];
        double vali = YGInti[];
        array_append (arrUnk, &vali, sizeof(double));
      }

#ifdef SOLVE_TEMPERATURE
      // Temperature field
      {
        double vali = TInt[];
        array_append (arrUnk, &vali, sizeof(double));
      }
#endif

      /**
      Gather parameters that must be passed to the
      non-linear system function. */

      UserDataNls data;

      coord o = {x,y,z};
      foreach_dimension()
        data.c.x = o.x;

      /**
      Solve the non-linear system of equations. */

      fsolve (EquationsGsl, arrUnk, &data);

      /**
      Recover Nls solution. */
      {
        double * unk = (double *)arrUnk->p;
        //double * unk = (double *) array_shrink (arrUnk);
        int count = 0;

        // mInt field
        for (int jj=0; jj<NLS; jj++) {
          scalar mEvapi = mEvapList[LSI[jj]];
          mEvapi[] = unk[count++];
        }

        // YLInt field
        for (int jj=0; jj<NLS; jj++) {
          scalar YLInti = YLIntList[jj];
          YLInti[] = unk[count++];
        }

        // YGInt field
        for (int jj=0; jj<NGS; jj++) {
          scalar YGInti = YGIntList[jj];
          YGInti[] = unk[count++];
        }
#ifdef SOLVE_TEMPERATURE
        // Temperature field
        {
          TInt[] = unk[count++];
        }
#endif
      }
      array_free (arrUnk);
    }
  }
}

#ifdef SOLVE_TEMPERATURE

#ifndef RADIATION_INTERFACE
# define RADIATION_INTERFACE 0.
#endif

double divq_rad_int (double TInti, double Tbulk = 300., double alphacorr = 1.) {
  //return alphacorr*5.669e-8*(pow(Tbulk, 4.) - pow(TInti, 4.));
  return alphacorr*5.670373e-8*(pow(Tbulk, 4.) - pow(TInti, 4.));
}

void EqTemperature (const double * xdata, double * fdata, void * params) {
  UserDataNls * data = (UserDataNls *)params;

  foreach_point (data->c.x, data->c.y, data->c.z, serial) {

    double TInti = xdata[0];
    bool success = false;

    double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true,  TInti, &success);
    double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TInti, &success);

#ifdef VARPROP
    double yl[NLS], yg[NGS];
    double xl[NLS], xg[NGS];
    double MWl[NLS], MWg[NGS];

    for (int jj=0; jj<NLS; jj++) {
      scalar YL = YLList[jj];
      yl[jj] = YL[];
      MWl[jj] = inMW[LSI[jj]];
    }
    mass2molefrac (xl, yl, MWl, NLS);

    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList[jj];
      yg[jj] = YG[];
      MWg[jj] = inMW[jj];
    }
    mass2molefrac (xg, yg, MWg, NGS);
#endif

    double vapheat = 0.;
    for (int jj=0; jj<NLS; jj++) {
      scalar mEvap = mEvapList[LSI[jj]];
#ifdef VARPROP
      scalar dhevjj = dhevList[jj];
      vapheat -= mEvap[]*dhevjj[];
#else
      vapheat -= mEvap[]*dhev;
#endif
    }

    fdata[0] = vapheat
        - divq_rad_int (TInti, TG0, RADIATION_INTERFACE)
#ifdef VARPROP
         + lambda1v[]*ltrgrad
         + lambda2v[]*gtrgrad
#else
         + lambda1*ltrgrad
         + lambda2*gtrgrad
#endif
         ;
  }
}

int EqTemperatureGsl (const gsl_vector * x, void * params, gsl_vector * f) {
 double * xdata = x->data;
  double * fdata = f->data;

  EqTemperature (xdata, fdata, params);
  return GSL_SUCCESS;
}


void ijc_CoupledTemperature ()
{
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      Array * arrUnk = array_new();
      {
        double vali = TInt[];
        array_append (arrUnk, &vali, sizeof(double));
      }
      UserDataNls data;

      coord o = {x,y,z};
      foreach_dimension()
        data.c.x = o.x;

#ifdef USE_GSL
      fsolve (EqTemperatureGsl, arrUnk, &data);
#endif
      {
        double * unk = (double *)arrUnk->p;
        TInt[] = unk[0];
        //if (unk[0] > 0. && unk[0] < 5000.) // The solution is reasonable
        //  TInt[] = unk[0];
      }
      array_free (arrUnk);
    }
  }
}


void EqBoilingTemperature (const double * xdata, double * fdata, void * params) {
  double Tb = xdata[0];
  double * xc = (double *)params;
  double sumKeq = 0.;
  for (int jj=0; jj<NLS; jj++) {
#ifdef USE_ANTOINE
    scalar YL = YLList[jj];
    sumKeq += YL.antoine (Tb, Pref)*xc[jj];
#endif
#ifdef USE_ANTOINE_OPENSMOKE
    sumKeq += opensmoke_antoine (Tb, Pref, jj)*xc[jj];
#endif
  }
  fdata[0] = sumKeq - 1.;
}

int EqBoilingTemperatureGsl (const gsl_vector * x, void * params, gsl_vector * f) {
  double * xdata = x->data;
  double * fdata = f->data;

  EqBoilingTemperature (xdata, fdata, params);
  return GSL_SUCCESS;
}

double boilingtemperature (double Tfg, double * xc)
{
  Array * arrUnk = array_new();
  array_append (arrUnk, &Tfg, sizeof(double));
#ifdef USE_GSL
  fsolve (EqBoilingTemperatureGsl, arrUnk, xc);
#endif
  double * unk = arrUnk->p;
  double Tb = unk[0];
  array_free (arrUnk);
  return Tb;
}

#endif // SOLVE_TEMPERATURE

