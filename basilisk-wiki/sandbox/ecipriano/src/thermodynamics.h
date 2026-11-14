/**
# Thermodynamics functions

Collection of thermodynamics functions, useful for
the evaporation models.

## Physical Constants
*/

#ifndef THERMODYNAMICS_H
# define THERMODYNAMICS_H

#define R_GAS 8.3144621 // Ideal gas constant [J/mol/K]

/**
## Clausius-Clapeyron
The function *clapeyron()* returns the thermodynamic
equilibrium constant: $P_{vap}(T)/P$ for a single chemical species.

### *clapeyron()*

* *Tint*: interface temperature [K]
* *Tboil*: boiling temperature [K]
* *dhev*: Enthalpy of evaporation [J/kg]
* *MWi*: Molecular weight [kg/kmol]
*/

double clapeyron (double Tint, double Tboil, double dhev, double MWi) {
  return exp (-dhev/R_GAS*MWi/1000. * (1./Tint - 1./Tboil));
}

/**
## Antoine

We define a function pointer for a generic antoine function.
The pointer must be set by the user, to a specific antoine
function depending on the chemical species under investigation.
*/

attribute {
  double (* antoine) (double, double);
}

/**
For consistency, the temperature and pressure values provided to
the antoine functions should be in SI units: K, Pa. Conversions
necessary for the specific Antoine equation are manages inside
each specific function. Each of the following functions return
the thermodynamic equilibrium constant $P_{vap}(T)/P$ using
the Antoine equation, whose parameters were taken from the
[NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C142825&Mask=4&Type=ANTOINE&Plot=on)
database.
*/

/**
### *antoine_heptane(T,P)*: Antoine equation for n-heptane.
*/

double antoine_heptane (double T, double P) {
  double A, B, C;
  if (T > 295)
    A = 4.02832, B = 1268.636, C = -56.199;
  else
    A = 4.81803, B = 1635.409, C = -27.338;
  return pow (10., A - B/(T + C)) / (P*1.e-5);
}

/**
### *antoine_decane(T,P)*: Antoine equation for decane.
*/

double antoine_decane (double T, double P) {
  double A, B, C;
  if (T > 367)
    A = 4.07857, B = 1501.268, C = -78.67;
  else
    A = 0.21021, B = 440.616, C = -156.896;
  return pow (10., A - B/(T + C)) / (P*1.e-5);
}

/**
### *antoine_hexadecane(T,P)*: Antoine equation for dodecane.
*/

double antoine_dodecane (double T, double P) {
  double A = 4.10549, B = 1625.928, C = -92.839;
  return pow (10., A - B/(T + C)) / (P*1.e-5);
}

/**
### *antoine_hexadecane(T,P)*: Antoine equation for hexadecane.
*/

double antoine_hexadecane (double T, double P) {
  double A = 4.17312, B = 1845.672, C = -117.054;
  return pow (10., A - B/(T + C)) / (P*1.e-5);
}

/**
### *antoine_methanol(T,P)*: Antoine equation for methanol.
*/

double antoine_methanol (double T, double P) {
  double A, B, C;
  if (T > 354)
    A = 5.15853, B = 1569.613, C = -34.846;
  else
    A = 5.20409, B = 1581.341, C = -33.50;
  return pow (10., A - B/(T + C)) / (P*1.e-5);
}

/**
## Composition Utilities

Functions that convert between mass and mole fractions, and that compute the
molecular weight of the mixture. They operate on plain arrays of fractions and
molecular weights, without any dependency on the thermodynamic state or on the
properties backend, therefore they belong here rather than in
[variable-properties.h](variable-properties.h).
*/

/**
## *mass2molefrac()*: Compute mole fractions from mass fractions

* *X*: vector filled with mole fractions
* *W*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

void mass2molefrac (double * X, const double * W, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += W[i]/MW[i];
  }
  for (int i=0; i<NS; i++) {
    X[i] = W[i]/MW[i]/(MWmix + 1.e-10);
  }
}

/**
## *mole2massfrac()*: Compute mass fractions from mole fractions

* *W*: vector filled with mole fractions
* *X*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

void mole2massfrac (double * W, const double * X, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += X[i]*MW[i];
  }
  for (int i=0; i<NS; i++) {
    W[i] = X[i]*MW[i]/(MWmix + 1.e-10);
  }
}

/**
## *mass2mw()*: Compute mixture molecular weight from mass fractions

* *W*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

double mass2mw (const double * W, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += W[i]/MW[i];
  }
  return 1./(MWmix + 1.e-10);
}

/**
## *mole2mw()*: Compute mixture molecular weight from mole fractions

* *X*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

double mole2mw (const double * X, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += X[i]*MW[i];
  }
  return MWmix;
}

/**
## *correctfrac()*: Close to 1 a vector of mass or mole fractions

* *X*: vector with mass or mole fractions
* *NS* total number of species (vector length)
*/

void correctfrac (double * X, const int NS)
{
  double sum = 0.;
  for (int i=0; i<NS; i++)
    sum += (X[i] >= 0.) ? X[i] : 0.;
  for (int i=0; i<NS; i++)
    X[i] = (X[i] >= 0.) ? X[i]/(sum + 1.e-10) : 0.;
}

#endif
