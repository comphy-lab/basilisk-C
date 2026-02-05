/**
# Temperature Profile Module Example
This is a simple example demonstrating how to use the Temperature Profile Module.
It sets up a temperature profile based on predefined time and temperature arrays,
and retrieves temperature values at specific times using linear interpolation.
*/

#include "temperature-profile.h"
double TG0 = 300.;

int main () {
  double timeprofile[] = {0, 1, 2, 3, 4, 5, 13.8996139, 41.6988417, 88.03088803, 166.7953668, 
    254.8262548, 342.8571429, 454.0540541, 574.5173745, 694.980695, 806.1776062, 
    917.3745174, 1037.837838, 1200};
  double temperatureprofile[] = {300, 309.492891, 366.8388626, 424.1848341, 486.7440758, 559.7298578, 
    611.8625592, 656.1753555, 697.8815166, 723.9478673, 736.9810427, 752.6208531, 
    750.014218, 752.6208531, 752.6208531, 750.014218, 750.014218, 750.014218, 750.014218};
  
  TemperatureProfile_Set(timeprofile, temperatureprofile, sizeof(timeprofile)/sizeof(double));

  TemperatureProfile_Print();

  fprintf(stderr, "T at t = 300: %g\n", TemperatureProfile_GetT(300));
  return 0;
}