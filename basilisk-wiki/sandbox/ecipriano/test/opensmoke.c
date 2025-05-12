/**
# OpenSMOKE Interface

We use the [opensmoke.h](../src/opensmoke.h) module to read
a kinetic scheme and to print out info, such as the number
of chemical species, their names, the molecular weights, ecc.).
The goal of this test is to check whether opensmoke is
correctly installed and compiled, all the relative paths
are correctly set and it works along with basilisk.
*/

#include "run.h"
#include "opensmoke.h"

int main (void) {

  /**
  We set the name of the kinetics folder path, relative to
  the directory $OPENSMOKE_INTERFACE/kinetics. We call the
  run function because [opensmoke.h](../src/opensmoke.h)
  automatically reads the kinetics in a defaults event.
  */

  kinfolder = "one-step/n-heptane";
  run();
}

/**
We print the information that we need to check. */

event logfile (i = 0) {
  fprintf (stderr, "Number of gas phase species = %d\n",
      OpenSMOKE_NumberOfSpecies());
  fprintf (stderr, "Number of liq phase species = %d\n",
      OpenSMOKE_NumberOfLiquidSpecies());
  fprintf (stderr, "Number of gas phase reactions = %d\n",
      OpenSMOKE_NumberOfReactions());

  for (int i=0; i<OpenSMOKE_NumberOfSpecies(); i++)
    fprintf (stderr, "Species[%d] = %s - molecular weight = %f [kg/kmol]\n",
        i, OpenSMOKE_NamesOfSpecies (i), OpenSMOKE_MW (i));
  fprintf (stderr, "Index of O2 is = %d\n",
      OpenSMOKE_IndexOfSpecies("O2"));

  for (int i=0; i<OpenSMOKE_NumberOfLiquidSpecies(); i++)
    fprintf (stderr, "Liquid Species[%d] = %s\n",
        i, OpenSMOKE_NamesOfLiquidSpecies (i));
}

