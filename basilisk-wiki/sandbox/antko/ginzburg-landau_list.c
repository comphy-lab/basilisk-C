/**
# The Ginzburg--Landau equation: a "full" implicitation

We solve the [complex Ginzburg--Landau
equation](/src/examples/ginzburg-landau.c) using [a solver for a list of equations](./list_solve.h). */
#include "grid/multigrid.h"
#define LIST_LEN 2 // The size of the system is defined here
#include "./list_solve.h"

/**
## Parameters

We use the same parameters as in the [Ginzburg-Landau example](/src/examples/ginzburg-landau.c). */
scalar Ar[], Ai[], A2[];

int main() {

  double beta = 1.5;
  size (100);
  init_grid (256);

  foreach() {
    Ar[] = 1e-4*noise();
    Ai[] = 1e-4*noise();
  }
  double dt = 5e-2;
  int i = 0;
  TOLERANCE = 1e-6;
/**
## Time integration 
The time loop is embedded in the main function. This is because the current implementation of the macros prevents the use of `run.h`.
*/
  for (double t = 0; t <= 150; t += dt, i++) {
    scalar r[], lambda[], rhsAr[], rhsAi[];
    foreach()
      A2[] = sq(Ar[]) + sq(Ai[]);
    foreach() {
      r[] = A2[]*beta;
      lambda[] = 1. - A2[];
      rhsAr[] = Ar[];
      rhsAi[] = Ai[];
    }
    list_solve (Ar, 
      (double[]){Ar[] - dt * (Ar[1] + Ar[-1] + Ar[0,1] + Ar[0,-1] - 4. * Ar[]) / sq(Delta) - dt * Ar[] * lambda[] - dt * r[] * Ai[], Ai[] - dt * (Ai[1] + Ai[-1] + Ai[0,1] + Ai[0,-1] - 4. * Ai[]) / sq(Delta) - dt * Ai[] * lambda[] + dt * r[] * Ar[]}, 
      (double[]){rhsAr[],rhsAi[]}, 
      (double[]){1. + 4. * dt / sq(Delta) - dt * lambda[], 1. + 4. * dt / sq(Delta) - dt * lambda[]}, 
      Ai);

  }
  FILE * fp2 = fopen ("dumpA2", "w");
  FILE * fpi = fopen ("dumpAi", "w");
  double lastx = HUGE;
  foreach() {
    if (x != lastx && lastx != HUGE)
    fputs ("\n", fp2); 
    fprintf (fp2, "%g %g %g\n", x, y, A2[]);
    lastx = x;
  }
  fclose (fp2);
  lastx = HUGE;
  foreach() {
    if (x != lastx && lastx != HUGE)
    fputs ("\n", fpi); 
    fprintf (fpi, "%g %g %g\n", x, y, Ai[]);
    lastx = x;
  }
  fclose (fpi);
}


/**
## Outputs

~~~gnuplot Amplitude (left) and phase (right) at t = 150
set term pngcairo size 600,600 transparent enhanced font "Helvetica,10"
set output 'sol.png'
unset key
unset colorbox
unset xtics
unset ytics
unset xlabel
unset ylabel
set view map
set size ratio -1
set pm3d map 
set palette gray
set multiplot layout 1, 2
splot 'dumpA2' using 1:2:3 w pm3d
splot 'dumpAi' using 1:2:3 w pm3d
unset multiplot
~~~
*/
