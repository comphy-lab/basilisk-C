/**
# Hydrostatic pressure 

Checks that we recover exact hydrostatic balance for two superposed
fluids with a density contrast of 100. 

~~~gnuplot Pressure
set xlabel 'z'
set ylabel 'Phi'
set xrange [0:0.3]
set yrange [0:0.9]
g = 9.81
rho1 = 1.
rho2 = 0.01
H0 = 0.29
hi = 0.087
plot 'log' w p t '', rho2*g*(H0 - x), rho2*g*(H0 - hi) + rho1*g*(hi - x)
~~~
*/

#include "grid/multigrid1D.h"
#define a_baro(eta, i) 0.
#include "layered/hydro.h"
#include "nh1.h"

double * rho = (double [60]){ 1. };

const double H0 = 0.29;
const double hi = 0.087;
double ai;

int main()
{
  size (6.);
  origin (- L0/2.);
  N = 8;
  nl = 60;
  G = 9.81;
  DT = 0.01;
#if 0
  lambda_b[] = {HUGE}; // free slip
#endif

  for (int l = 0; l < nl; l++)
    rho[l] = l < nl/2 ? 1. : 0.01;
  
  ai = 0.0;
  run();
}

/**
The initial layer positions. */

event init (i = 0)
{
  foreach() {
    double zi = ai*x/(L0/2.);
    foreach_layer()
      h[] = point.l < nl/2 ? (hi - zi)/(nl/2) : (H0 - hi + zi)/(nl/2);
  }
}

event logfile (i = 1; i <= 10)
{
  if (i == 0)
    fprintf (stderr, "\n\n# ai = %g\n", ai);
  foreach() {
    double z = zb[];
    foreach_layer() {
      // the pressure is defined at the bottom of the layer
      fprintf (stderr, "%g %g\n", z, phi[]);
      z += h[];
    }
  }
}
