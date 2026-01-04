/**
# Cantera

Test for the cantera interface. We check and provide an usage example of:

* How to read a kinetics in yaml format, processed from chemkin using
  [ck2yaml](https://cantera.org/3.1/userguide/ck2yaml-tutorial.html)
* How to compute variable material properties
* How to integrates stiff ideal-gas constant pressure batch reactors

These operations are limited to the gas phase.
*/

#include "common.h"
#include "run.h"

/**
We include the main interface with cantera, which is responsible for reading the
kinetics, the variable material properties, and the chemical reactions module.
*/

#include "cantera/cantera.h"
#include "cantera/properties.h"
#include "cantera/chemistry.h"

/**
We need to store the number of species and the species names. */

int NS, NLS;
char ** species = NULL, ** species_liq = NULL;

/**
#fixme: liquid properties should not be necessary. */

double Dmix1 = 0.;
double cp1 = 1000.;
double dhev = 2050.;
double lambda1 = 0.3;
double mu1 = 1e-3;
double rho1 = 1000.;

/**
We read the kinetics, store the number of species and their names, run all the
events, and free the heap allocations. */

int main (void) {
  kinetics ("cantera.yaml", &NS);
  kinetics_liquid ("cantera.yaml", &NLS);

  fprintf (stderr, "NS = %d\n", NS);
  fprintf (stderr, "NLS = %d\n", NLS);

  species = new_species_names (NS);
  species_liq = new_species_names_liquid (NLS);

  run();

  free_species_names (NS, species);
  free_species_names (NLS, species_liq);
  kinetics_clean();
  return 0;
}

event init (i = 0) {
  /**
  We print the molecular weight and the mixture properties. */

  double MW[NS];
  molecular_weights (NS, MW);
  for (int i = 0; i < NS; i++)
    fprintf (stderr, "MW[%d] = %g\n", i, MW[i]);

  for (int i = 0; i < NS; i++)
    fprintf (stderr, "species[%d] = %s\n", i, species[i]);

  for (int i = 0; i < NLS; i++)
    fprintf (stderr, "species_liquid[%d] = %s\n", i, species_liq[i]);

  double x[NS];
  for (int i = 0; i < NS; i++)
    x[i] = 0.;
  x[index_species ("CH3OH")] = 0.4;
  x[index_species ("O2")] = 0.126;
  x[index_species ("N2")] = 0.474;

  ThermoState tsg;
  tsg.T = 1200., tsg.P = 101325., tsg.x = x;

  fprintf (stderr, "gas density = %g\n", tp2.rhov (&tsg));
  fprintf (stderr, "gas viscosity = %g\n", tp2.muv (&tsg));
  fprintf (stderr, "gas conductivity = %g\n", tp2.lambdav (&tsg));
  fprintf (stderr, "gas heat capacity = %g\n", tp2.cpv (&tsg));

  double D[NS], cp[NS], betaT[NS];
  tp2.diff (&tsg, D);
  tp2.cps (&tsg, cp);
  tp2.betaT (&tsg, betaT);

  for (int i = 0; i < NS; i++)
    fprintf (stderr, "D[%s] = %g\n", species[i], D[i]);

  for (int i = 0; i < NS; i++)
    fprintf (stderr, "cp[%s] = %g\n", species[i], cp[i]);

  /**
  We integrate a batch reactor with combustion chemistry. */

  double y0[NS+1];
  for (int i = 0; i < NS; i++)
    y0[i] = x[i];
  y0[NS] = tsg.T;

  double s0[NS+1];
  UserDataODE data;
  data.T = tsg.T;
  data.P = tsg.P;
  data.sources = s0;

  double t = 0.;
  while (t < 1e-3) {
    fprintf (stderr, "batch %g %g", t, y0[NS]);
    for (int i = 0; i < NS; i++)
      fprintf (stderr, " %g", y0[i]);
    fprintf (stderr, "\n");

    stiff_ode_solver (NULL, NS+1, 1e-5, y0, &data);
    t += 1e-5;
  }
}

/**
## Results

~~~gnuplot Temperature profiles
set grid
set xlabel "time [s]"
set ylabel "Temperature [K]"
plot "<grep batch log" u 2:3 w l t "Temperature"
~~~

~~~gnuplot Species profiles
set xlabel "time [s]"
set ylabel "Mass Fractions [-]"
plot "<grep batch log" u 2:4 w l t "CH3OH", \
     "<grep batch log" u 2:5 w l t "O2", \
     "<grep batch log" u 2:6 w l t "CO2", \
     "<grep batch log" u 2:7 w l t "CO", \
     "<grep batch log" u 2:8 w l t "H2O", \
     "<grep batch log" u 2:9 w l t "N2"
~~~
*/
