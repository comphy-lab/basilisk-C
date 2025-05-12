/**
# Boiling Temperature

we test the calculation of the boiling temperature for a binary
mixture of n-heptane/n-hexadecane at ambient pressure and at
different composition. */

/**
We set data for the multicomponent model. */

#define NGS 3
#define NLS 2

#define USE_GSL 0
#define USE_ANTOINE
#define SOLVE_TEMPERATURE

#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "multicomponent.h"

char* gas_species[NGS] = {"NC7H16", "NC16H34", "N2"};
char* liq_species[NLS] = {"NC7H16", "NC16H34"};
char* inert_species[1] = {"N2"};
double gas_start[NGS] = {0.0, 0.0, 1.0};
double liq_start[NLS] = {0.5, 0.5};
double inDmix1[NLS] = {1e-6, 1e-6};
double inDmix2[NGS] = {1e-5, 1e-5, 1e-5};
double inKeq[NLS] = {0.8, 0.4};
double lambda1, lambda2, cp1, cp2, dhev, TG0, TL0;

int main (void) {
  run();
}

event init (i = 0) {

  /**
  We set the `antoine` functions to the liquid mass fraction
  fields. */

  scalar C7  = YLList[0];
  scalar C16 = YLList[1];

  C7.antoine = antoine_heptane;
  C16.antoine = antoine_hexadecane;

  /**
  We loop over 5 different compositions and we call the function that
  computes the boiling temperature using 300 K as first guess value.
  */

  double comp[] = {0., 0.25, 0.5, 0.75, 1.};
  for (int i=0; i<5; i++) {
    double xc[] = {comp[i], 1.-comp[i]};
    fprintf (stderr, "x = %.2f - Tb [K] = %.2f\n", comp[i], boilingtemperature (300., xc));
  }
}

/**
As expected, the boiling temperature decreases with increasing
amount of n-heptane in the mixture:

~~~
x = 0.00 - Tb [K] = 559.94
x = 0.25 - Tb [K] = 426.52
x = 0.50 - Tb [K] = 397.01
x = 0.75 - Tb [K] = 381.67
x = 1.00 - Tb [K] = 371.58
~~~
*/
