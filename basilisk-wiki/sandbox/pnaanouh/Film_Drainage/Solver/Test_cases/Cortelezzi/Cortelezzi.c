/** 
# Viscous capillary waves - Comparison with Prosperetti & Cortelezzi solution.

We test here the ability of the multilayer solver and the 
surface tension extension to reproduce correctly viscous capillary waves 
for various regimes.
*/

/** 
## Includes
We use a constant-resolution grid with the multilayer the Navier--Stokes solver 
and the [surface tension additive](../../hydro-tension.h). */
#include "grid/multigrid1D.h"
#include "../hydro-tension.h"
#include "../nh.h"
#include "../viscous_surface.h"
#include "layered/remap.h"

/**
##Parameters */
#define Bo 1.
#define Oh 1.


#define lay 2
#define LEVEL 7

#define delta 0.05 [1]

double T0 = 5.E5 [0,1];
double k = 0.01[-1];

/**
## Initialisation */
int main() {
  
  
  N = 1 << LEVEL;
  nl = lay;
  L0 = 2.*pi/k;

  nu = Oh*1 [2,-1];
  G = Bo*1 [1,-2];
  
  TOLERANCE = 1e-8 [*];
  CFL_H = 0.1;
  h_diffusion=true;
  
  run();
 
}

/**
The initial condition is a small amplitude plane wave of wavelength
unity and uniform surface tension and viscosity. */

event init (t = 0) {
  periodic(right);
  foreach () {
    double H = 1. [1] + delta*cos(k*x);
    foreach_layer () {
      h[] = H/nl;
    }
  }
}

/**
##Outputs */

void plot_profile (double t, FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "p [0:%g][-0.04:1.1]'-' u 1:3:2 w filledcu lc 3 t ''\n", t,L0);
  foreach()
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n\n");
  fflush (fp);
}
event animate_profile (i += 100.)
{
  
static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
if (i == 0)
   fprintf (fp, "set term x11\n");
plot_profile (t, fp);
  
}


/**
We output the evolution of the film (amplitude and vorticity)
 in a file indexed with the simulation parameters. */
event amplitude (t += T0/2000) {
  double max = statsf(eta).max;
  char name[80];
  if (visc_activate)
    sprintf (name, "ampSimu_Oh_%g_k_%g_cor_N_%d_L_%d", nu, k, N, nl);
  else
    sprintf (name, "ampSimu_Oh_%g_k_%g_N_%d_L_%d", nu, k, N, nl);
  static FILE * fpw = fopen (name, "w");
  fprintf (fpw, "%g %g\n", t, (max-1.)/delta);
  fflush (fpw);
}

event end(t=T0)
{}
