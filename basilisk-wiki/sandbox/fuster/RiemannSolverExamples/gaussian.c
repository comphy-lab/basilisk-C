/**
# Shock tube with Riemann solver

We solve the Euler equations for a compressible gas. We also need to
compute volume fractions for the initial condition. */

#include "grid/multigrid.h"
#include "compressible.h"

double tend = 3.;

int main() {

  /**
  We make boundary conditions free outflow. */

  foreach_dimension() {
    w.n[right] = neumann(0);
    w.n[left]  = neumann(0);
  }
  
  size(20);
  X0 = -L0/2.;
  N = 128;
  run(); 
}

event init (t = 0)
{ 
  
  double cson = sqrt(gammao);

  foreach() {
    double p = 1.+ 1.e-3*exp(-x*x);
    rho[] = 1. + (p - 1.)/sq(cson);
    foreach_dimension()
      w.x[] = 0.;
    E[] = p/(gammao - 1.);
  }
}

event endprint (t = tend) {
  
  foreach () {
    double xref = fabs(x) - tend*sqrt(gammao);
    printf("%g %g \n", xref, (E[]*(gammao-1)-1.)/1.e-3);
  }
  printf("\n");

}

/**

~~~gnuplot Pressure profile
set xrange[-3:3]
p "out" u 1:2 t 'Numerical' w lp pt 7, 0.5*exp(-x*x) t 'Theory' w l lw 2 lc 0
~~~ 

*/
