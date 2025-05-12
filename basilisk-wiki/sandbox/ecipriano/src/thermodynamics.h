/**
# Thermodynamics functions

Collection of thermodynamics functions, useful for
the evaporation models.

## Physical Constants
*/

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

#ifdef USE_ANTOINE

attribute {
  double (* antoine) (double, double);
}

#endif

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

