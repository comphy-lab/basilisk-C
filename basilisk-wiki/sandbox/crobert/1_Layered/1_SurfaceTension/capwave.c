/**
# Capillary wave 

This is an adaptation of the classical test case first proposed 
in [Popinet & Zaleski, 1999](/src/references.bib#popinet1999).*/

/** 
## Includes
We use a constant-resolution grid with the multilayer the Navier--Stokes solver 
and the [surface tension additive](/../tension.h). */
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"
#include "capwave_Cortelezzi.h"

/**
We will store the accumulated error in *se* and the number of samples
in *ne*. */
double se = 0; int ne = 0;

/**
## Initialisation */
#define k 0.5
int main() {
  /**
  The surface tension is one and the viscosity is constant. */
  L0 = 2.*pi/k;
  periodic(right);
  
  G = 0;
  nu = 0.005;
  corr_dux = true;
  corr_dwx = true;

  /**
  We vary the resolution to check for convergence. */
  nl = 32;
  for (N = 4; N <= 64; N *= 2) {
    se = ne = 0;
    run();
  }
  
  N = 64;
  for (nl = 1; nl <= 32; nl *= 2) {
    se = ne = 0;
    run();
  }
}

/**
The initial condition is a small amplitude plane wave of wavelength
unity. */
event init (t = 0) {
  foreach () {
    double H = 1. + 0.01*cos (k*x);
    foreach_layer () {
      h[] = H/nl;
      sigma[] = 1.;
    }
  }
}

/**
We output the amplitude at times matching exactly those in the
reference file. */
event amplitude (t += 0.130931; t <= 78.4278) {
  double max = (statsf(eta).max - 1.)/0.01;

  /**
  We output the corresponding evolution in a file indexed with the
  number of grid points *N*. */

  char name[80];
  sprintf (name, "wave-N%d-nl%d", N, nl);
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t, max);
  fflush (fp);

  /**
  To compute the RMS error, we get data from the reference file
  *prosperetti.h* and add the difference to the accumulated error. */

  se += sq(max - cortelezzi[ne][1]); ne++;
}

/**
At the end of the simulation, we output on standard error the
resolution (number of grid points per wavelength) and the relative RMS
error. */

event error (t = end)
  fprintf (stderr, "%g %d %g\n", N/L0, nl, sqrt(se/ne));

/**
## Results

~~~gnuplot Evolution of the amplitude of the capillary wave as a function of non-dimensional time $\tau=\omega_0 t$
set xlabel 'tau'
set ylabel 'Relative amplitude'
plot '../capwave_Cortelezzi.h' u 2:4 w l t "Cortelezzi", \
     'wave-N32-nl32' every 10 w p t "Basilisk"
~~~

~~~gnuplot Convergence of the RMS error as a function of resolution (number of grid points per wavelength)
set xlabel 'Number of grid points'
set ylabel 'Relative RMS error'
set logscale y
set logscale x 2
plot [0.3:5]'log' every 1::0::4 u 1:3 t "Multilayer" w lp, 1./x t "First order"
~~~

~~~gnuplot Convergence of the RMS error as a function of resolution (number of layers)
set xlabel 'Number of layers'
set ylabel 'Relative RMS error'
set logscale y
set logscale x 2
plot [2:50]'log' every 1::5::10 u 2:3 t "Multilayer" w lp, 1./x t "First order"
~~~
*/