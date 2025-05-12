/**
This example aims to provide documentation about the beahviour of the dirichlet and the radiation boundary conditions when of the form of a temporal sinus.
*/

/**
~~~gnuplot Dirichlet conditions
reset session
set term gif animate delay 5
set output 'animation_0.gif'
unset key
getValue (row, col, filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')
data_out = "data_out_0"
set xr[getValue (1, 2, data_out):getValue (1, 2, data_out)+getValue (2, 2, data_out)] # X0 and L0 are, respectively, at rows 1 and 2, columns 2 of "data_out"
set yr[-1:1]
ftime = 'time_0'
do for [i=1:getValue (7, 2, data_out)] { # counter is at row 7, column 2 of "data_out"
  set title getValue (i, 1, ftime)
  file = sprintf ('beach_0_%05.0f', i)
  plot file using 1:2:3 with filledcurves lc rgb '#56B4E9', \
  file using 1:3 with filledcurves above y1=getValue (8, 2, data_out) lc rgb '#000000', \
  file using 1:2 with lines lc rgb '#0000FF', \
  file using 1:4 with lines lc rgb '#55DD55'
}
~~~
*/

/**
~~~gnuplot Radiation conditions
reset session
set term gif animate delay 5
set output 'animation_1.gif'
unset key
getValue (row, col, filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')
data_out = "data_out_1"
set xr[getValue (1, 2, data_out):getValue (1, 2, data_out)+getValue (2, 2, data_out)] # X0 and L0 are, respectively, at rows 1 and 2, columns 2 of "data_out"
set yr[-1:1]
ftime = 'time_1'
do for [i=1:getValue (7, 2, data_out)] { # counter is at row 7, column 2 of "data_out"
  set title getValue (i, 1, ftime)
  file = sprintf ('beach_1_%05.0f', i)
  plot file using 1:2:3 with filledcurves lc rgb '#56B4E9', \
  file using 1:3 with filledcurves above y1=getValue (8, 2, data_out) lc rgb '#000000', \
  file using 1:2 with lines lc rgb '#0000FF', \
  file using 1:4 with lines lc rgb '#55DD55'
}
~~~
*/

/**
~~~gnuplot Custom radiation conditions
reset session
set term gif animate delay 5
set output 'animation_2.gif'
unset key
getValue (row, col, filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')
data_out = "data_out_2"
set xr[getValue (1, 2, data_out):getValue (1, 2, data_out)+getValue (2, 2, data_out)] # X0 and L0 are, respectively, at rows 1 and 2, columns 2 of "data_out"
set yr[-1:1]
ftime = 'time_2'
do for [i=1:getValue (7, 2, data_out)] { # counter is at row 7, column 2 of "data_out"
  set title getValue (i, 1, ftime)
  file = sprintf ('beach_2_%05.0f', i)
  plot file using 1:2:3 with filledcurves lc rgb '#56B4E9', \
  file using 1:3 with filledcurves above y1=getValue (8, 2, data_out) lc rgb '#000000', \
  file using 1:2 with lines lc rgb '#0000FF', \
  file using 1:4 with lines lc rgb '#55DD55'
}
~~~
*/

/**
As we can see, in the long run radiation conditions cause the water level to go higher than with the Dirichlet conditions.
*/

/**
# Code
*/


/**
We include the same scripts as in `beach.c`, but here the value of `ML` gives the number of layers.
*/
#define ML 3

#include "grid/multigrid1D.h"
#if ML
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#else
#include "green-naghdi.h"
#endif



#define NUMBER_OF_TESTS 3


/**
We define the final time `Tfinal` and the level of refinement ($N=2^{LEVEL}$).
*/
#define Tfinal 450.

#define LEVEL 10



/**
We can control wether we want an animation to be displayed while the script is running and/or an animation to be saved after the script is done, with gnuplot.\
The value is the number of frames between snapshots.
*/
#define animation_while_running 0
#define animation 25
#if animation
int counter = 0;
FILE * ftime;
#endif



/**
We keep some parameters form `beach.c` :
the amplitude of the wave `A`, the nominal water depth `h0` and the `slope`.
*/
double A = 0.3;
double h0 = 1.;
double slope = 1./20.;

scalar H[];


/**
We define the domain as `[X0;-X0+20]`.
We get the domain from `beach.c` with `X0=-40` aproximately.
*/

int test = 0; // 0 = Dirichlet; 1 = Radiation; 2 = Custom radiation; 3 = Custom Lande radiation
int main ()
{
  X0 = -40.;
  L0 = -X0 + 20.;
  N = 1 << LEVEL;
  G = 1.;
  
  #if ML
  nl = ML;
  breaking = 0.07;
  #else
  alpha_d = 1.;
  breaking = 0.4;
  #endif

  DT = 1.;
  
  int test_number;
  for (test_number = 0; test_number < NUMBER_OF_TESTS; test_number++) {
    test = test_number;
    #if animation
    counter = 0;
    #endif
    run();
  }
}



/**
## Initialization
*/
event init (i = 0)
{
  double c = sqrt(G*(1. + A)*h0);
  foreach() {
    zb[] = max (-h0, slope*x);
    #if ML
    double eta = 0.01*A*(x < X0 + 5. ? sin(pi*(x - X0)/5.) : 0.); // Small sinus wave because if not the multilayer solver behaves bizarrely
    foreach_layer() {
      h[] = max (0., eta - zb[])/(double)nl;
      u.x[] = c*eta/(h0 + eta);
    }
    #else
    double eta = 0.;
    h[] = max (0., eta - zb[]);
    u.x[] = c*eta/(h0 + eta);
    #endif
  }

  /**
  Boundary conditions, generating simple waves.\
  We test Dirichlet, radiation, and custom radiation, the latter being the radiation function but without taking into account the actual water depth but the at rest water depth in place.
  */
  if (test == 0)
    u.n[left] = dirichlet (A*sin(t*0.5));
  else if (test == 1)
    u.n[left] = radiation (A*sin(t*0.5));
  else if (test == 2)
    u.n[left] = sqrt(G*(max(0., 1.))) - sqrt(G*(max(0., A*sin(t*0.5) - zb[0]))); // Custom "radiation"
}


/**
## Events
Friction event from `beach.c`.\
As this event is done at every step, we take this oportunity to compute the depth field `H`, used in other places in the code.
This will also be useful for clarity, because we won't have to write `#if ML...#else...#endif` every time.
*/
event friction (i++) {
  foreach() {
    #if ML
    double Q = 0.;
    H[] = 0.;
    foreach_layer() {
      H[] += h[];
      Q += h[]*u.x[];
    }
    double a = H[] < dry ? HUGE : 1. + 5e-3*dt*fabs(Q)/sq(H[]);
    foreach_layer()
      u.x[] /= a;
    #else
    H[] = h[];
    double a = h[] < dry ? HUGE : 1. + 5e-3*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
    #endif
  }
}




/**
Event from `beach.c` to display an animation of the wave while the programm is running.
*/
#if animation_while_running
event running_animation (i += animation_while_running) {
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "p [%g:%g][-1.:1.]'-' u 1:3:2 w filledcu lc 3 t '',"
	   " '' u 1:(-1):3 t '' w filledcu lc -1\n", t, X0, X0+L0);
  foreach(serial)    
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n\n");
}
#endif



/**
If we want an animation at the end, we save the quantities we want periodically.
*/
#if animation
//event animation_datafile (i += animation) {
event animation_datafile (t = 350.; t += 1.) {
  char filename[30];
  counter += 1;
  sprintf(filename, "beach_%1.1i_%5.5d", test, counter);
  FILE * fanim = fopen (filename, "w");
  foreach(serial) {
    #if ML
    fprintf (fanim, "%g %g %g %g\n", x, eta[], zb[], H[]);
    #else
    fprintf (fanim, "%g %g %g %g\n", x, eta[], zb[], h[]);
    #endif
  }
  fclose (fanim);
  if (counter == 1) {
    char filename_time[7];
    sprintf (filename_time, "time_%1.1i", test);
    ftime = fopen (filename_time, "w");
  }
  fprintf (ftime, "%g\n", t);
}
#endif






/**
In the final event, we write a file, called `data_out`, to save useful information about the implementation, particularly for the animation script.
*/

event end (t = Tfinal) {
  #if animation
  fclose (ftime);
  #endif
  char filename[11];
  sprintf (filename, "data_out_%1.1i", test);
  FILE * fp = fopen (filename, "w");
  // "Physical parameters"
  fprintf (fp, "X0= %g\n", X0);
  fprintf (fp, "L0= %g\n", L0);
  // Numerical parameters
  fprintf (fp, "N= %i\n", N);
  fprintf (fp, "nl= %i\n", ML);
  // End state
  fprintf (fp, "i= %i\n", i);
  fprintf (fp, "t= %g\n", t);
  // Animation data
  #if animation
  fprintf (fp, "counter= %i\n", counter);
  fprintf (fp, "min(zb)= %g\n", statsf(zb).min);
  #endif
  fclose (fp);  
}
