/**
# Optically Thin Radiation Model

We test the calculation of the radiation using
the optically-thin model at different temperautres.

We define variables that are necessary for the
radiation model to compile.
*/

#define NGS 2
const char* gas_species[2] = {"CO2", "H2O"};

#include "common.h"
#include "utils.h"
#include "radiation.h"

int main (void) {

  /**
  The properties for the radiation model
  are gathered in the following struct. */

  OpticallyThinProperties otp;
  otp.T = 1400.;
  otp.P = 101325.;
  otp.xH2O = 0.21;
  otp.xCO2 = 0.1;

  /**
  We set the radiation model, which is
  set to no_radiation by default. */

  divq_rad = optically_thin;

  /**
  We loop over different temperatures and
  compute the $\nabla\cdot\mathbf{q}_{rad}$
  contribution. */

  for (otp.T = 300.; otp.T <= 3000.; otp.T += 50.) {
    fprintf (stderr, "divq_rad = %f\n", divq_rad(&otp));
  }
}
