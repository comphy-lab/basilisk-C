/**
This code is a simulation of the Plateau-Rayleigh instability. 

The simulation does not stop at MAXTIME. It's stop before due to an arithmetic
error. This error occur in the library geometry.h, arround the line 40, in the
function line_alpha. It is due to the fact that we want to compute the square
root of a negative number.

There is 2 way of avoiding this bug.

The first one lie in the liquid parameter. Indeed, if we use the air-water
parameter, with a small value of Oh, it works.

The second one is to disable the adaptativity of the mesh.

The bug occurs fastly enough: it takes only one minute for the code to reach it.

We think that this bug is due to some incompatibility between surface-tension
and adaptativity. */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#define RHOR 0.966
#define MUR 100.

#define LEVELINI 6
#define REFINEMENTMAX 9
#define REFINEMENTMIN 5
#define MAXTIME 60.
#define L0 10

double R = 1, epsilon = 0.1;
double Oh = 0.560;

int main(){
  f.sigma = 1.;
  mu1 = Oh;
  mu2 = mu1/MUR;
  rho1 = 1.;
  rho2 = rho1/RHOR;
  size (L0);
  init_grid (1 << LEVELINI);
  run();
}


/** We setup a cylinder of liquid with a perturbation on it, to initiate the
break-up process. */ 

event init (t = 0) {
  fraction (f, R*(1. + epsilon*cos((M_PI*x/10))) - y);
}

/**
If we disable the adaptation in the code, there is no more arithmetic exception.
*/
event adapt (i++; t<=MAXTIME) {
  adapt_wavelet ({f, u}, (double[]){0.001,0.2,0.2,0.2},
  REFINEMENTMAX, REFINEMENTMIN);
  fprintf (stderr, "%d %g %d %d %d %ld\n", i, t, 
  mgp.i, mgpf.i, mgu.i, grid->tn);
}