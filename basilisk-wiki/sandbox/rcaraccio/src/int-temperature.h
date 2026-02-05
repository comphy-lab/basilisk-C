/**
# Interface temperature solver
This header contains the implementation of the solver for the interface temperature.
We exploit the GSL library to solve the non-linear equation for the interface temperature.
*/

#include "intgrad.h"
#ifndef USE_GSL
# define USE_GSL 1
#endif
#include "fsolve-gsl.h"

#ifndef RADIATION_TEMP
#define RADIATION_TEMP TG0
#endif

//Extern variables
extern scalar fS, fG;
extern face vector fsS, fsG;
extern vector lambda1v, lambda2v;
extern double TG0;
extern scalar TInt, TS, TG;

const double char_emissivity = 0.1;
const double wood_emissivity = 0.6;

typedef struct {
  coord c;
} UserDataNls;

#ifndef RADIATION_INTERFACE
# define RADIATION_INTERFACE 0.
#endif

// Radiative heat flux at the interface
double divq_rad_int (double TInti, double Tbulk = 300., double alphacorr = 1.) {
  return alphacorr*5.670373e-8*(pow(Tbulk, 4.) - pow(TInti, 4.));
}

// Non-linear equation for the interface temperature to be zeroed
int EqTemperature (const gsl_vector * xdata, void * params, gsl_vector * fdata) {
  UserDataNls * data = (UserDataNls *)params;
  double TInti = gsl_vector_get(xdata, 0);
  
  Point point = locate(data->c.x, data->c.y, data->c.z);
  // foreach_point(data->c.x, data->c.y, data->c.z, serial) {

  bool success = false;

  double gradTGn = ebmgrad(point, TG, fS, fG, fsS, fsG, true, TInti, &success);
  double gradTSn = ebmgrad(point, TS, fS, fG, fsS, fsG, false, TInti, &success);

  coord n = facet_normal(point, fS, fsS);
  normalize(&n);
  n.x = fabs(n.x); n.y = fabs(n.y);

  double lambda1vh = n.x / (n.x + n.y) * lambda1v.x[] + n.y / (n.x + n.y) * lambda1v.y[];
  double lambda2vh = n.x / (n.x + n.y) * lambda2v.x[] + n.y / (n.x + n.y) * lambda2v.y[];

  scalar YCHAR = YSList[OpenSMOKE_IndexOfSolidSpecies("CHAR")];
  double char_fraction = YCHAR[] / fS[];
  double emissivity = char_fraction * char_emissivity + (1. - char_fraction) * wood_emissivity;

  gsl_vector_set(fdata, 0,
                 -divq_rad_int(TInti, RADIATION_TEMP, emissivity) 
                 + lambda1vh * gradTSn 
                 + lambda2vh * gradTGn);
  // }
  return GSL_SUCCESS;
}

void ijc_CoupledTemperature() {

  foreach() {
    if (f[]>F_ERR && f[] < 1.-F_ERR) {
      gsl_vector *unk = gsl_vector_alloc(1);
      gsl_vector_set(unk, 0, TInt[]);

      UserDataNls data;
      coord o = {x,y,z};
      foreach_dimension()
        data.c.x = o.x;

      fsolve_gsl (EqTemperature, unk, &data, "EqTemperature");

      TInt[] = gsl_vector_get(unk, 0);
      gsl_vector_free(unk);
    }
  }
}
