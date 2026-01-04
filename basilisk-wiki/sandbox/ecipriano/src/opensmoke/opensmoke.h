/**
# OpenSMOKE++

Basilisk interface for the [OpenSMOKE++ library](#cuoci2015opensmoke++),
which is a C++ framework for numerical simulations of reacting systems
with detailed kinetic mechanisms, including thousands of chemical
species and reactions. It provides the functionalities for managing
the calculation of thermodynamics and transport properties of mixtures
and for handling the chemical kinetics.

## OpenSMOKEInterface
This module requires an external library, which contains both the
OpenSMOKE++ library itself, and the C interface which contains the
functions that can be called from Basilisk. It can be downloaded
from [github](https://github.com/edocipriano/OpenSMOKEppInterface)
and it can be easily compiled using cmake, provided the following
external dependencies: 
[OpenSMOKE++](#cuoci2015opensmoke++), 
[boost](https://www.boost.org/), 
[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), 
[libconfig++](https://github.com/hyperrealm/libconfig).
*/

#ifndef OPENSMOKE_H
#define OPENSMOKE_H
#define OPENSMOKE 1

#include "OpenSMOKE_Interface.h"

#pragma autolink -L$OPENSMOKE_INTERFACE/build -lopensmoke

/**
## User Data

The variable *kinfolder*, with the path of the kinetic scheme folder
relative to the path $OPENSMOKE_INTERFACE/kinetics, must be set by the
user.
*/

char * kinfolder, * liqfolder;

void molecular_weights (size_t ns, double * MW) {
  for (int i = 0; i < ns; i++)
    MW[i] = OpenSMOKE_MW (i);
}

int index_species (char * name) {
  return OpenSMOKE_IndexOfSpecies (name);
}

char ** new_species_names (size_t ns) {
  char ** species = (char **)malloc (ns*sizeof (char *));
  for (size_t i = 0; i < ns; i++)
    species[i] = strdup (OpenSMOKE_NamesOfSpecies (i));
  return species;
}

char ** new_species_names_liquid (size_t ns) {
  char ** species = (char **)malloc (ns*sizeof (char *));
  for (int i = 0; i < ns; i++) {
    const char * name = OpenSMOKE_NamesOfLiquidSpecies (i);
    int len = strlen (name);
    char corrname[len+1];
    strcpy (corrname, name);
    corrname[3 <= len ? len-3 : 0] = '\0';
    species[i] = strdup (corrname);
  }
  return species;
}

void free_species_names (size_t ns, char ** species) {
  for (size_t i = 0; i < ns; i++)
    free (species[i]);
  free (species);
}

#define MAX_KINFOLDER_LEN 1200

void kinetics (char * kinfolder, int * NS = NULL) {
  char kinfolder_root[MAX_KINFOLDER_LEN];
  sprintf (kinfolder_root, "%s/kinetics/%s/kinetics",
      getenv ("OPENSMOKE_INTERFACE"), kinfolder);

  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder_root);

  if (NS) *NS = OpenSMOKE_NumberOfSpecies();
}

void kinetics_liquid (char * kinfolder, int * NS = NULL) {
  char kinfolder_root[MAX_KINFOLDER_LEN];
  sprintf (kinfolder_root, "%s/kinetics/%s/kinetics",
      getenv ("OPENSMOKE_INTERFACE"), kinfolder);

  OpenSMOKE_ReadLiquidKinetics (kinfolder_root);

  if (NS) *NS = OpenSMOKE_NumberOfLiquidSpecies();
}

void properties_liquid (char * liqfolder) {
  char liqfolder_root[MAX_KINFOLDER_LEN];
  sprintf (liqfolder_root, "%s/kinetics/LiquidProperties/%s",
      getenv ("OPENSMOKE_INTERFACE"), liqfolder);

  OpenSMOKE_ReadLiquidProperties (liqfolder_root);
}

void kinetics_clean (void) {
  OpenSMOKE_Clean();
}

#endif

/**
## References

~~~bib
@article{cuoci2015opensmoke++,
  title={OpenSMOKE++: An object-oriented framework for the numerical modeling of reactive systems with detailed kinetic mechanisms},
  author={Cuoci, Alberto and Frassoldati, Alessio and Faravelli, Tiziano and Ranzi, Eliseo},
  journal={Computer Physics Communications},
  volume={192},
  pages={237--264},
  year={2015},
  publisher={Elsevier}
}
~~~
*/

