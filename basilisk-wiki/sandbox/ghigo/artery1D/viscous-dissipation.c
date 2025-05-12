/**
# Viscous dissipation

We solve the propagation of a wave using the viscous 1D blood flow equations in a straight artery. As the wave propagates, viscous dissipation modifies the amplitude of the wave. */

/**
# Asymptotic viscous dissipation

For a small amplitude $QIN$ of the flow rate signal, an asymptotic solution can be obtained, where the envelop of the viscous dissipation is described by the following expression:

$$
QIN * \exp\left(-\frac{CF}{2}x \right)
$$
*/

/**
# Code
*/

#include "grid/cartesian1D.h"
#include "bloodflow.h"

double R0 = 1.;
double K0 = 1.e4;
double L = 200.;

double CF = 0.1;

double T0 = 1.;
double SH = 1.e-3;
double AIN = 0., QIN = 0.;

double celerity (double a, double k) {

  return sqrt(0.5*k*sqrt(a));
}

double inlet (double t, double t0) {

  if (t < t0) return max(0., 0.5*(1. + cos(pi + 2.*pi/T0*t)));
  else return 0.;
}

int main() {

  origin (0., 0.);
  L0 = L;

  AIN = pi*pow(R0*(1 + SH), 2.);
  QIN = SH*AIN*celerity(AIN, K0);

  N = 2048;
    
  run();
  
  return 0; 
}

q[left] = dirichlet (QIN*inlet (t, T0));

event defaults (i=0) {

  gradient = order1;
  riemann = hll_glu;
}

event init (i=0) {
  
  foreach() {
    k[] = K0;
    a0[] = pi*pow(R0, 2.);
    a[] = a0[];
    q[] = 0.;
  }
}

event viscosity (i++, last) {

  scalar sa[], sq[];
  foreach() {
    sa[] = 0.;
    sq[] = -CF*q[]/a[];
    q[] += dt*sq[];
  }
  boundary({sa,sq});
}

event field (t = {0.5, 1., 1.5, 2., 2.5}) {

  foreach() {
    fprintf (stderr, "%g, %.6f, %.6f\n", x, a[], q[]) ;
  }
  fprintf (stderr, "\n\n") ;

}

event end (t = 2.5) {
  printf ("#Viscous dissipation test case completed\n") ;
}

/**
# Plots
*/

/**
~~~gnuplot Spatial evolution of the cross-sectional area $A$ computed at $t={0.5, 1, 1.5, 2, 2.5}$ for $N=128$ 

reset

mycolors = "dark-blue red sea-green dark-violet orange"
mypoints = '1 2 3 4 6'

set datafile separator ','
set style line 1 lt -1 lw 2 lc 'black'

set output 'A.png'
set xlabel 'x [cm]'
set ylabel 'A [cm^2]'
set key center center

plot 'log' i 0 u 1:2 every 5 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 2 lw 1 t 't=0.5', \
     'log' i 1 u 1:2 every 5 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=1', \
     'log' i 2 u 1:2 every 5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 2 lw 1 t 't=1.5', \
     'log' i 3 u 1:2 every 5 w p lt word(mypoints, 4) lc rgb word(mycolors, 4) ps 2 lw 1 t 't=2', \
     'log' i 4 u 1:2 every 5 w p lt word(mypoints, 5) lc rgb word(mycolors, 5) ps 2 lw 1 t 't=2.5'
~~~

~~~gnuplot Spatial evolution of the flow rate $Q$ computed at $t={0.5, 1, 1.5, 2, 2.5}$ for $N=128$ 

k0 = 1.e4
a0 = pi
l = 200.
cf = 0.1
sh = 1.e-3
ain = a0*(1 + sh)**2.;
qin = sh*ain*sqrt(0.5*k0*sqrt(ain));
amp(x) = qin*exp(-cf/2.*x/l)

set output 'Q.png'
set xlabel 'x [cm]'
set ylabel 'Q [cm^3/s]'
set key right center

plot amp(x) w l ls 1 t 'Asymptotic', \
     'log' i 0 u 1:3 every 5 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=1', \
     'log' i 1 u 1:3 every 5 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=1', \
     'log' i 2 u 1:3 every 5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 2 lw 1 t 't=1.5', \
     'log' i 3 u 1:3 every 5 w p lt word(mypoints, 4) lc rgb word(mycolors, 4) ps 2 lw 1 t 't=2', \
     'log' i 4 u 1:3 every 5 w p lt word(mypoints, 5) lc rgb word(mycolors, 5) ps 2 lw 1 t 't=2.5'
~~~
*/
