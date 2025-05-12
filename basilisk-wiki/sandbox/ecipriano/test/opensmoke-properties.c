
//#include "run.h"
//#include "utils.h"
#include "navier-stokes/centered.h"
#include "opensmoke-properties.h"
#include "two-phase.h"

int main (void) {
  mu1 = 1. , mu2 = 1.;
  kinfolder = "skeletal/methanol";
  run();
}

/**
## Properties

We update the material properties according to the
thermodynamic state of the mixture.
*/

event properties (i++) {

  /**
  We first set the thermodynamic state of the mixture: */

  double x[OpenSMOKE_NumberOfSpecies()];
  double q[OpenSMOKE_NumberOfLiquidSpecies()];
  for (int jj=0; jj<OpenSMOKE_NumberOfSpecies(); jj++)
    x[jj] = 0.;
  x[OpenSMOKE_IndexOfSpecies("O2")] = 0.21;
  x[OpenSMOKE_IndexOfSpecies("N2")] = 0.79;
  q[0] = 1.;

  ts1.T = 300., ts1.P = 101325., ts1.x = q;
  ts2.T = 300., ts2.P = 101325., ts2.x = x;

  /**
  We update the properties of the variable fields. */

  fprintf (stderr, "Gas Properties\n");
  fprintf (stderr, "rho2 = %f\n", tp2.rhov (&ts2));
  fprintf (stderr, "mu2 = %f\n", tp2.muv (&ts2));
  fprintf (stderr, "lambda2 = %f\n", tp2.lambdav (&ts2));
  fprintf (stderr, "cp2 = %f\n", tp2.cpv (&ts2));

  for (unsigned int i=0; i<OpenSMOKE_NumberOfSpecies(); i++) {
    fprintf (stderr, "Dmix[%s] = %f\n",
        OpenSMOKE_NamesOfSpecies (i), tp2.diff (&ts2, i));
  }

  fprintf (stderr, "\n");

  fprintf (stderr, "Liquid Properties\n");
  fprintf (stderr, "rho1 = %f\n", tp1.rhov (&ts1));
  fprintf (stderr, "mu1 = %f\n", tp1.muv (&ts1));
  fprintf (stderr, "lambda1 = %f\n", tp1.lambdav (&ts1));
  fprintf (stderr, "cp1 = %f\n", tp1.cpv (&ts1));

  for (unsigned int i=0; i<OpenSMOKE_NumberOfLiquidSpecies(); i++) {
    fprintf (stderr, "pvap[%s] = %f\n",
        OpenSMOKE_NamesOfLiquidSpecies (i), tp1.pvap (&ts1, i));
    fprintf (stderr, "dhev[%s] = %f\n",
        OpenSMOKE_NamesOfLiquidSpecies (i), tp1.dhev (&ts1, i));
    fprintf (stderr, "sigma[%s] = %f\n",
        OpenSMOKE_NamesOfLiquidSpecies (i), tp1.sigmas (&ts1, i));
  }
}

event logfile (i++) {
  for (double Temp = 290.; Temp <= 500.; Temp += 5.) {
    double x[] = {1};
    ts1.T = Temp, ts1.P = 101325., ts1.x = x;
    fprintf (stdout, "%g %g\n", Temp, liqprop_thermal_expansion (&tp1, &ts1));
  }
}

