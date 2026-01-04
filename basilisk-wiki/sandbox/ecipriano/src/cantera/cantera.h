/**
# Cantera

Basilisk interface for [Cantera](https://cantera.org/index.html), which
which is a C++ framework for numerical simulations of reacting systems with
detailed kinetic mechanisms, including thousands of chemical species and
reactions. It provides the functionalities for managing the calculation of
thermodynamics and transport properties of mixtures and for handling the
chemical kinetics.

This module requires the [installation](https://cantera.org/stable/install/index.html)
of the Cantera library.
*/

#ifndef CANTERA_H
#define CANTERA_H
#define CANTERA 1

#include <sys/stat.h>
#include <cantera/clib/ct.h>
#include <cantera/clib/ctreactor.h>
#include <cantera/clib/ctfunc.h>
#include <cantera/clib/ctmultiphase.h>
#include <cantera/clib/ctonedim.h>
#include <cantera/clib/ctrpath.h>
#include <cantera/clib/ctsurf.h>

#define CANTERA 1

#pragma autolink -lcantera -lstdc++ -lfmt -lyaml-cpp -lblas -llapack -lsundials_arkode -lsundials_cvode -lsundials_cvodes -lsundials_ida -lsundials_idas


int soln = -1;
int thermo = -1, kin = -1, tran = -1;
int soln_liq = -1;
int thermo_liq = -1, kin_liq = -1, tran_liq = -1;

#define MAX_KINFOLDER_LEN 1200

void kinetics (char * kinfolder, int * NS = NULL) {
  struct stat sb;
  if (stat (kinfolder, &sb) == 0)
    soln = soln_newSolution (kinfolder, "gas", "default");
  else {
    char kinfolder_root[MAX_KINFOLDER_LEN];
    sprintf (kinfolder_root, "%s/kinetics/%s/kinetics/kinetics.yaml",
        getenv ("OPENSMOKE_INTERFACE"), kinfolder);
    soln = soln_newSolution (kinfolder_root, "gas", "default");
  }

  if (soln < 0) {
    fprintf (stderr, "error: Cantera was unable to read the kinetics.\n");
    exit(1);
  }
  thermo = soln_thermo (soln);
  kin = soln_kinetics (soln);
  tran = soln_transport (soln);

  if (NS) *NS = thermo_nSpecies (thermo);
}

void kinetics_liquid (char * kinfolder, int * NS = NULL) {
  struct stat sb;
  if (stat (kinfolder, &sb) == 0) {
    soln_liq = soln_newSolution (kinfolder, "liq", "default");
  }
  else {
    char kinfolder_root[MAX_KINFOLDER_LEN];
    sprintf (kinfolder_root, "%s/kinetics/%s/kinetics/kinetics.yaml",
        getenv ("OPENSMOKE_INTERFACE"), kinfolder);

    soln_liq = soln_newSolution (kinfolder_root, "liq", "default");
  }
  if (soln < 0) {
    fprintf (stderr, "error: Cantera was unable to read the liquid kinetics.\n");
    exit(1);
  }
  thermo_liq = soln_thermo (soln_liq);
  kin_liq = soln_kinetics (soln_liq);
  tran_liq = soln_transport (soln_liq);

  if (NS) *NS = thermo_nSpecies (thermo_liq);
}

void properties_liquid (char * liqfolder) {
  return;
}

void kinetics_clean (void) {
  ct_appdelete();
}

void molecular_weights (size_t ns, double * MW) {
  thermo_getMolecularWeights (thermo, ns, MW);
}

int index_species (char * name) {
  return thermo_speciesIndex (thermo, name);
}

char ** new_species_names (size_t ns) {
  char ** species = (char **)malloc (ns*sizeof (char *));
  for (size_t i = 0; i < ns; i++) {
    species[i] = (char *)malloc (80*sizeof (char));
    thermo_getSpeciesName (thermo, i, 80, species[i]);
  }
  return species;
}

char ** new_species_names_liquid (size_t ns) {
  char ** species = (char **)malloc (ns*sizeof (char *));
  char fullspecies[80];
  for (size_t i = 0; i < ns; i++) {
    species[i] = (char *)malloc (80*sizeof (char));
    thermo_getSpeciesName (thermo_liq, i, 80, fullspecies);
    size_t len = strlen (fullspecies);
    memcpy (species[i], fullspecies, len-3);
    species[i][len-3] = '\0';
  }
  return species;
}

void free_species_names (size_t ns, char ** species) {
  for (size_t i = 0; i < ns; i++)
    free (species[i]);
  free (species);
}

#endif
