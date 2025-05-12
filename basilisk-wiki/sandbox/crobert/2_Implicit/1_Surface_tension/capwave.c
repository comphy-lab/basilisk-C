/**
# Capillary wave 

This is an adaptation of the classical test case first proposed 
in [Popinet & Zaleski, 1999](/src/references.bib#popinet1999).*/

/** 
## Includes
We use a constant-resolution grid with the multilayer the Navier--Stokes solver 
and the [surface tension additive](/crober/2_Implicit/hydro-tension.h). */
#include "grid/multigrid1D.h"
#include "crobert/2_Implicit/hydro-tension.h"
#include "layered/nh.h"
#include "crobert/2_Implicit/viscous_surface.h"
#include "layered/remap.h"
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
  The surface tension is one by default and the viscosity is constant. */
  L0 = 2.*pi/k;
  periodic(right);
  TOLERANCE = 1e-8;
  linearised = true;
  CFL_H = 0.5;
  
  G = 0;
  nu = 0.005;

  /**
  We vary the resolution to check for convergence. */
  nl = 32;
  for (N = 4; N <= 64; N *= 2) {
    se = ne = 0;
    run();
  }
  
  N = 64;
  for (nl = 1; nl <= 16; nl *= 2) {
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
  To compute the relative RMS error, we get data from the reference file
  *capwave_Cortelezzi.h* and add the difference to the accumulated error.
  The reference was computed thanks to the Python code 
  [cortelezzi_analytical.py](cortelezzi_analytical.py) inspired from the work
  of [Cortelezzi & Propseretti](#prosperetti1982)*/

  se += sq(max - cortelezzi[ne][1]); ne++;
}

/**
At the end of the simulation, we output on standard error the
resolution (number of grid points per wavelength) and the relative RMS
error. */

event error (t = end)
  fprintf (stderr, "%d %d %g\n", N, nl, sqrt(se/ne));

/**
## Results

~~~gnuplot Evolution of the amplitude of the capillary wave as a function of non-dimensional time $\tau=\omega_0 t$
set xlabel 'tau'
set ylabel 'Relative amplitude'
plot '../capwave_Cortelezzi.h' u 2:4 w l t "Cortelezzi", \
     'wave-N64-nl32' every 10 w p t "Basilisk"
~~~

~~~gnuplot Convergence of the RMS error as a function of resolution (number of grid points)
set xlabel 'Number of grid points'
set ylabel 'Relative RMS error'
set logscale y
set logscale x 2
plot 'log' every 1::0::4 u 1:3 t "Multilayer" w lp, 1./x t "First order"
~~~

~~~gnuplot Convergence of the RMS error as a function of resolution (number of layers)
set xlabel 'Number of layers'
plot [1:20]'log' every 1::5::9 u 2:3 t "Multilayer" w lp, 1./x t "First order"
~~~

Accuracy is first-order in N but saturates for high N, even with a low tolerance.

## References
~~~bib
@article{prosperetti1982,
  title={Small-amplitude waves produced by a submerged vorticity distribution on the surface of a viscous liquid},
  author={Prosperetti, Andrea and Cortelezzi, L},
  journal={The Physics of Fluids},
  volume={25},
  number={12},
  pages={2188--2192},
  year={1982},
  publisher={American Institute of Physics}
}
~~~
*/
