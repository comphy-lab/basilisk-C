/**
# Radiation Models

Radiation models for combustion simulations. Currently,
just the optically-thin model is implemented.
This implementation is based on the information reported
in [https://tnfworkshop.org/radiation/](https://tnfworkshop.org/radiation/)
*/

#define RADIATION

/**
We define the Stefan-Boltzmann constant $[W/m^2/K^4]$. */

#define STEFAN_BOLTZMANN  5.669e-8

/**
## Optically-Thin Model

This is the simplest model to implement, it computes the
heat loss for radiation according to the following equation:
$$
\nabla\cdot \mathbf{q}_{rad} =
-4\sigma\sum\limits_{i=1}^{NGS}a_{p,i} P x_i (T^4 - T_b^4)
$$
where $\sigma$ is the Stefan-Boltzmann constant, $a_{p,i}$
is the Planck mean absorption coefficient of species $i$,
$P x_i$ is the partial pressure of the chemical species,
$T$ is the local flame temperature, while $T_b$ is a baseline
temperature (approximated as 300K).

The following properties are required by the radiation
model:

* *T*: temperature of the system [K]
* *P*: pressure of the system [Pa]
* *xCO2*: mole fraction of CO2 [moli/molt]
* *xH2O*: mole fraction of H2O [moli/molt]

Only CO2 and H2O are considered, since their contribution to
the radiation is the gratest, if compared with the other
species.
*/

typedef struct {
  double T, P;
  double xCO2, xH2O;
  double xCO, xCH4;
} OpticallyThinProperties;

/**
We also declare a struct that gathers optically-thin settings. */

struct OpticalyThinModel {
  int indexCO2, indexH2O;
  int indexCO, indexCH4;
};

struct OpticalyThinModel otm = {
  .indexCO2 = -1,
  .indexH2O = -1,
  .indexCO  = -1,
  .indexCH4 = -1,
};

/**
We compute the radiation according to the optically-thin model. */

double optically_thin (void * p) {

  /**
  We extract the data from the void pointer. We limit the temperature
  to 2500K because it's the range of reliability of the Plank-mean
  absorption coefficients. Normally the temperature is not higher than
  that, but when the spark is present it can give problems. */

  OpticallyThinProperties * otp;
  otp = (OpticallyThinProperties *)p;

  double T = otp->T, P = otp->P, uT = 1000./T;
  double xCO2 = otp->xCO2, xH2O = otp->xH2O;
  double xCO = otp->xCO, xCH4 = otp->xCH4;

  /**
  The pressure value in SI (Pa) is converted in atm. */

  //P *= 9.86923e-6;

  /**
  We compute the Plank-mean absorption using coefficients
  from ([Chu 2014](#chu2014calculations)), that do not
  degenerate above 2500K. All the coefficients have units
  [1/m/bar]. */

  double apCO2 = 18.741 +uT*(-121.31+uT*(273.5 +uT*(-194.05 +uT*( 56.31 + uT*(-5.8169)))));
  double apH2O = -0.23093 + uT*(-1.1239 + uT*(9.4153 + uT*(-2.9988 + uT*(0.51382 + uT*(-1.8684e-5)))));
  double apCO  = (T < 750.) ? (4.7869 + T*(-0.06953 + T*(2.95775e-4 + T*(-4.25732e-7 + T*2.02894e-10)))) :
    (10.09 + T*(-0.01183 + T*(4.7753e-6 + T*(-5.87209e-10 + T*-2.5334e-14))));
  double apCH4 = 6.6334 + T*(-0.0035686 + T*(1.6682e-08 + T*(2.5611e-10 - 2.6558e-14*T)));

  double sum_pa = (apCO2*xCO2 + apH2O*xH2O + apCO*xCO + apCH4*xCH4)*P/1e+5;

  /**
  We return the flux of radiant energy. */

  return -4.*STEFAN_BOLTZMANN*sum_pa*( pow (T, 4.) - pow (300., 4.) );
}

/**
## No Radiation

Implementation of a dummy no_radiation model, which returns
a 0 divq contribution from the radiation.
*/

double no_radiation (void * p) {
  return 0;
}

/**
## Main radiation function

We declare a function pointer to a generic radiation
model that computes the term $\nabla\cdot q_{rad}$
for the energy equation.
*/

double (* divq_rad) (void * p) = no_radiation;

/**
## Init Event

We overload the init event in order to setup properties
for the radiation models at the beginning of the simulation.
*/

event init (i = 0) {
  if (divq_rad == optically_thin) {

    for (int jj=0; jj<NGS; jj++) {
      if (strcmp (gas_species[jj], "CO2") == 0)
        otm.indexCO2 = jj;
      if (strcmp (gas_species[jj], "H2O") == 0)
        otm.indexH2O = jj;
      if (strcmp (gas_species[jj], "CO") == 0)
        otm.indexCO = jj;
      if (strcmp (gas_species[jj], "CH4") == 0)
        otm.indexCH4 = jj;
    }

    // Species not found
    if (otm.indexCO2 == -1)
      fprintf (ferr, "WARNING: Optically thin model did not found CO2\n");
    if (otm.indexH2O == -1)
      fprintf (ferr, "WARNING: Optically thin model did not found H2O\n");
    if (otm.indexCO == -1)
      fprintf (ferr, "WARNING: Optically thin model did not found CO\n");
    if (otm.indexCH4 == -1)
      fprintf (ferr, "WARNING: Optically thin model did not found CH4\n");
  }
}

/**
## References

~~~bib
@article{grosshandler1993radcal,
  title={Radcal: A narrow band model for radiation},
  author={Grosshandler, William L},
  journal={Calculations in a Combustion Environment, NIST Technical Note},
  volume={1402},
  year={1993}
}
@article{chu2014calculations,
  title={Calculations of narrow-band transimissities and the Planck mean absorption coefficients of real gases using line-by-line and statistical narrow-band models},
  author={Chu, Huaqiang and Gu, Mingyan and Zhou, Huaichun and Liu, Fengshan},
  journal={Frontiers in Energy},
  volume={8},
  pages={41--48},
  year={2014},
  publisher={Springer}
}
~~~

*/

