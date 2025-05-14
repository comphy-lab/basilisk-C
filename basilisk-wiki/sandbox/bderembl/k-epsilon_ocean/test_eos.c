/**
   Test the non linear eos-80 equation of state for sea water.


~~~gnuplot Sea Water Equation of State (EOS-80).
unset key
set size ratio 0.5
unset surface
set pm3d map
set contour base
set cntrparam levels 10
set cntrlabel onecolor
set cntrlabel format "%4.0f"
set cntrlabel start 5 interval 50
set xlabel "Salinity (PSU)"
set ylabel "Temperature (degC)"
splot 'out' u 3:4:5 w l lc rgbcolor "black", 'out' u 3:4:5 w labels boxed
~~~

 */

#define EOS_UNESCO 1
#include "eos_ocean.h"
#include "utils.h"
#include "output.h"

/**
   Bounds to compute the T-S diagram
 */

double T0 = 0;
double T1 = 30.0;

double S0 = 20;
double S1 = 40;

scalar T[];
scalar S[];
scalar rho[];
scalar b[];

int main(){

  N = 64;
  init_grid(N);

  foreach(){
    S[] = S0 + x*(S1-S0); 
    T[] = T0 + y*(T1-T0); 
    rho[] = eos_rho(T[],S[]);
  }

  output_field ({S, T, rho});


}
