#define LB_HP_FILTER 0

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension-ssvm.h"

/**
We make sure that the boundary conditions for the face-centered
velocity field are consistent with the centered velocity field (this
affects the advection term). */

uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

/**
We will store the accumulated error in *se* and the number of samples
in *ne*. */

double se = 0; int ne = 0;

int main() {

  /**
  The domain is 2x2 to minimise finite-size effects. The surface
  tension is one and the viscosity is constant. */
  
  double ratios[5] = {15., 8., 4., 2., 1.};
  ratio_dt_st = 15.;

  size (2. [1]);
  Y0 = -L0/2.;
  f.sigma = 1.;
  TOLERANCE = 1e-6 [*];
  mu1 = mu2 = 0.0182571749236;
  
  N = 128;
  se = 0, ne = 0;
  
  


  /**
  We vary the resolution to check for convergence. */
  
  for (int i = 0; i <= 4; i++) {
    ratio_dt_st = ratios[i];
    run();
  }
}


/**
The initial condition is a small amplitude plane wave of wavelength
unity. */

FILE * fp = NULL;

event init (t = 0) {

  char name[80];
  sprintf (name, "capwave-%d-%d", 7, (int) ratio_dt_st);
  if (fp)
    fclose (fp);
  fp = fopen (name, "w");
  
  double k = 2., a = 0.01;
  fraction (f, y - a*cos (k*pi*x));
}

/**
By default tracers are defined at $t-\Delta t/2$. We use the *first*
keyword to move VOF advection before the *amplitude* output i.e. at
$t+\Delta/2$. This improves the results. */

event vof (i++, first);


event amplitude (i++) {

  /**
  To get an accurate amplitude, we reconstruct interface position
  (using height functions) and take the corresponding maximum. */

  scalar pos[];
  position (f, pos, {0,1 [0]});
  double max = statsf(pos).max;

  /**
  We output the corresponding evolution in a file indexed with the
  number of grid points *N*. */

  fprintf (fp, "%g %g\n", t*11.1366559937, max);
  fflush (fp);
}

event end (t = 2.2426211256);


/**
## Results

~~~gnuplot Evolution of the amplitude of the capillary wave as a function of non-dimensional time $\tau=\omega_0 t$
set xlabel 'tau'
set ylabel 'Relative amplitude'
plot '../prosperetti.h' u 2:4 w l t "Prosperetti", \
     'capwave-7-1' w l t "{/Symbol D} t = 1 {/Symbol D} t_{ST}, CSVM", \
     'capwave-7-2' w l t "{/Symbol D} t = 2 {/Symbol D} t_{ST}, CSVM", \
     'capwave-7-4' w l t "{/Symbol D} t = 4 {/Symbol D} t_{ST}, CSVM", \
     'capwave-7-8' w l t "{/Symbol D} t = 8 {/Symbol D} t_{ST}, CSVM", \
     'capwave-7-15' w l t "{/Symbol D} t = 15 {/Symbol D} t_{ST}, CSVM"
     
~~~

*/