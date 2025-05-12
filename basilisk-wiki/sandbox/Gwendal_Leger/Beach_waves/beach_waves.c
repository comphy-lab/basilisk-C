/**
This is a small example heavily based on [beach.c](http://basilisk.fr/src/test/beach.c).\
It is basically a fusion of my [first](http://basilisk.fr/sandbox/Gwendal_Leger/wave_topo/beach_topo.c) 
and [second](http://basilisk.fr/sandbox/Gwendal_Leger/waves_1D/waves_1D.c) tests, in which one can specify the bathymetry *and* the type of waves used.
*/


/**
~~~gnuplot Sinus wave on a custom bathymetry
reset session
set term gif animate delay 5
set output 'animation.gif'
unset key
getValue (row, col, filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')
data_out = "data_out"
set xr[getValue (1, 2, data_out):getValue (1, 2, data_out)+getValue (2, 2, data_out)] # X0 and L0 are, respectively, at rows 1 and 2, columns 2 of "data_out"
set yr[-1:1]
ftime = 'time'
file = 'animation_data'
if (getValue (7, 2, data_out) != 0) {
do for [i=1:getValue (7, 2, data_out)] { # counter is at row 7, column 2 of "data_out"
  set title getValue (i, 1, ftime)
  plot \
  file using 1:2:3 index i with filledcurves lc rgb '#56B4E9', \
  file using 1:3 index i with filledcurves above y1=getValue (8, 2, data_out) lc rgb '#000000', \
  file using 1:2 index i with lines lc rgb '#0000FF'
}
} else {
set title "Dummy plot"
plot x title "Dummy plot"
}
~~~
*/

/**
We include the same scripts as in `beach.c`, but here the value of `ML` gives the number of layers.\
If we choose not to use the multilayer solver, we can choose between using Green-Naghdi or Saint-Venant.
*/
#define ML 2

#include "grid/multigrid1D.h"
#if ML
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#else
#define GN 0
#if GN
#include "green-naghdi.h"
#else
#include "saint-venant.h"
#endif
#endif


/**
`interpolation.h` is a custom script made by me.
It interpolates values in different manners, see the file directly for more information
and input the desired data in the `interpolation_data.h` file.\
`output_gnuplot1D.h` is a custom script with a routine to take a snapshot of the fluid at a specific moment.
*/
#include "interpolation.h"
#include "output_gnuplot1D.h"
#include "DFT.h"
#include "waves.h"


/**
We define the final time `Tfinal` and the level of refinement ($N=2^{LEVEL}$).
*/
#define Tfinal 150.

#define LEVEL 10


/**
We can control wether we want an animation to be displayed while the script is running and/or an animation to be saved after the script is done, with gnuplot.\
The value is the number of frames between snapshots.
*/
#define animation_while_running 0
#define animation 7
int counter = 0;
#if animation
FILE *ftime, *fanim;
#endif


/**
We can choose among the same waves than in [my first example](http://basilisk.fr/sandbox/Gwendal_Leger/wave_topo/beach_topo.c),
or by supplying a 0 we can use the wave theory :

* 1 = square
* 2 = sinus
* 3 = pointy hat
* 4 = pointy sinus
* 5 = full sinus
* 6 = parabola
* 7 = C^oo_o
* <0 = soliton (the original from the `beach.c` file)
* else = wave theory

*/
#define wave_choice 0 // Type of wave (1 = square; 2 = sinus; 3 = pointy hat; 4 = pointy sinus; 5 = full sinus; 6 = parabola; 7 = C^oo_o; <0 = soliton)


/**
We keep some parameters form `beach.c` :
the amplitude of the wave `A`, the nominal water depth `h0` and the `slope`.
*/
double A = 0.3;
double h0 = 1.;
double slope = 1./20.;

scalar H[], U[];
#if ML
scalar Zl; // Altitude of (the middle of) each layer
#else
scalar Zl[];
#endif


/**
We take a snapshot after some time, at time `t1`.
*/
const double t1 = 10.;


wave_spectrum waves;

/**
We define the domain as `[X0;-X0+20]`.
We get the domain from `beach.c` with `X0=-40` aproximately.
*/
int main()
{
  X0 = -40.;
  L0 = -X0 + 20.;
  N = 1 << LEVEL;
  G = 1.;
  
  DT = 0.5; // In case of implicit scheme >>> u=0 which leads to dt>>1
  
  #if ML
  nl = ML;
  breaking = 0.07;
  #else
  #if GN
  alpha_d = 1.;
  breaking = 0.4;
  #endif
  #endif
  /**
     If we chose to use it, we setup the wave spectrum that will be used, by using the data from `waves_input.h` file.
  */
  #if wave_choice == 0
  waves = init_waves ();
  #endif
  
  run();
}


/**
# Wave choices
We keep the functions necessary for the soliton wave from `beach.c`.
*/
double sech2 (double x) {
  return sq (2./(exp(x) + exp(-x)));
}
double soliton (double x, double t) {
  double c = sqrt(G*(1. + A)*h0), psi = x - c*t;
  double k = sqrt(3.*A*h0)/(2.*h0*sqrt(h0*(1. + A)));
  return A*h0*sech2 (k*psi);
}
/**
We also implement some "homemade" waves, all taking the space of the five leftmost units of the domain.
*/
double initial_wave_height (double x, double t) {
  /** "Square" wave. */
#if wave_choice == 1
  return A*(x < X0 + 5. ? 1. : 0.);
  /** Sinus wave. */
#elif wave_choice == 2
  return A*(x < X0 + 5. ? sin(pi*(x - X0)/5.) : 0.);
  /** "Pointy hat" wave. */
#elif wave_choice == 3
  return A*(x < X0 + 5. ? x < X0 + 2.5 ? (x - X0)/2.5 : -(x - X0 - 5.)/2.5 : 0.);
  /** "Pointy sinus" wave. */
#elif wave_choice == 4
  return A*(x < X0 + 5. ? (x < X0 + 1.25 ? -(x - X0)/1.25 : (x > X0 + 3.75 ? -(x - X0 - 5.)/1.25 : (x - X0 - 2.5)/1.25)) : 0.);
  /** "Full" sinus wave. */
#elif wave_choice == 5
  return A*(x < X0 + 5. ? sin(pi + pi*(x - X0)/2.5) : 0.);
  /** Parabola. */
#elif wave_choice == 6
  return A*(x < X0 + 5. ? -1./6.25*pow(x - X0, 2) + 5./6.25*(x - X0) : 0.);
  /** The classical example of a function that is infinetly differentiable with a compact support. */
#elif wave_choice == 7
  return A*(x < X0 + 5. ? exp(-1./(1. - pow(fabs(x - X0 - 2.5)/2.5, 2))) / exp(-1.) : 0.);
  /** The soliton from `beach.c`. */
#elif wave_choice < 0
  double k = sqrt(3.*A/4/cube(h0)), L = 2./k*acosh(sqrt(1./0.05));
  return soliton (x + h0/slope + L/2., t);
  /** Lake at rest. */
#else
  return 0.;
#endif
}


/**
# Initialization
*/
event init (i = 0)
{
  #if ML
  Zl = new scalar[nl];
  #endif
  /**
  Function from `interpolation.h` that reads a file and interpolate the data read to fill the bathymetry `zb`.
  The data is read in the `interpolation_data.h` file.
  */
  interpolate_bathymetry_from_datafile (zb);
  
  double c = sqrt(G*(1. + A)*h0);
  foreach() {
    eta[] = initial_wave_height (x, t);
    H[] = max (0., eta[] - zb[]);
    U[] = c*eta[]/(h0 + eta[]);
    double zl = zb[];
    #if ML
    foreach_layer() {
      h[] = H[]/(double)nl;
      u.x[] = U[];
      zl += h[]*0.5;
      Zl[] = zl;
      zl += h[]*0.5;
    }
    #else
    h[] = H[];
    u.x[] = U[];
    #endif
  }
   /**
  Because of the homogenous Neumann conditions, the free surface elevation will be too high with the square wave.
  This is due to the fact that the "tail" of the wave is not at level h=h0 but at h=h0+A.
  So we manually fix the water level of the first 2 cells to h=h0.
  */
  #if wave_choice == 1
  h[0] = h0; h[1] = h0;
  #endif
  
  /**
     We take a snapshot of the initial state by using the function from `output_gnuplot1D.h`.
  */
#define ZRANGE_SNAPSHOT -0.5, 0.5
  snapshot_fluid ("eta0", eta, zb, 0., ZRANGE_SNAPSHOT);

  /**
  Boundary conditions (for now, only with lake at rest initialization).
  */
  #if wave_choice == 0
  u.n[left] = VeloX (x, 0., Zl[] - zb[], t, waves, -zb[]);
  #else
  u.n[left] = neumann (0.);
  #endif
}


/**
# Events
*/

/**
We compute the depth field `H` and the (vertically) average velocity field `U` at every step, especially with the mutlilayer solver in mind.
This will also be useful for clarity, because we won't have to write `#if ML...#else...#endif` every time.
*/
event compute_vertical_scalars (i++) {
  foreach() {
    #if ML
    H[] = 0.;
    U[] = 0.;
    foreach_layer() {
      H[] += h[];
      U[] += u.x[]*h[];
    }
    U[] = H[] > dry ? U[]/H[] : 0.; // If the cell is dry, we prevent a division by zero by manually setting the cell velocity to zero
    #else
    H[] = h[];
    U[] = u.x[];
    #endif
  }
}


/**
Friction event from `beach.c`.
*/
event friction (i++) {
  foreach() {
    #if ML
    double Q = H[]*U[];
    double a = H[] < dry ? HUGE : 1. + 5e-3*dt*fabs(Q)/sq(H[]);
    U[] = 0.; // We need to re-compute the average velocity field
    foreach_layer() {
      u.x[] /= a;
      U[] += u.x[]*h[];
    }
    U[] = H[] > dry ? U[]/H[] : 0.;
    #else
    double a = h[] < dry ? HUGE : 1. + 5e-3*dt*norm(u)/h[];
    u.x[] /= a;
    U[] = u.x[];
    #endif
  }
}




/**
We save/display the iteration, time, time step, maximum and mean depth, and maximum and mean horizontal velocity.
*/
event logfile (i += 5) {
  double maxh = 0., meanh = 0., meanu = 0., maxu = 0.;
  maxh = statsf(H).max;
  meanh = statsf(H).sum/(double)N;
  maxu = statsf(U).max;
  meanu = statsf(U).sum/(double)N;
  fprintf (stderr, "%i %6.6g %3.3e  %3.3e %3.3e %3.3e %3.3e\n",
	   i, t, dt, maxh, meanh, maxu, meanu);
}






/**
We use Basilisk's tide gauges to measure the evolution of the free surface elevation during the simulation, and then perform a discrete Fourier transform of the signal.
 */
const int n_gauges = 3; // Number of tide gauges
Gauge gauges_custom[] = {
			  {"Gauge_left", 0., 0., "Left"},
			  {"Gauge_middle", 0., 0., "Middle"},
			  {"Gauge_right", 0., 0., "Right"},
			  {NULL}};
event init_gauges (i = 0) {
  int ii;
  for (ii = 0; ii < n_gauges; ii++)
    gauges_custom[ii].x = X0 + L0*((double)ii/(double)(n_gauges - 1) - (ii == n_gauges-1 ? 0.001 : 0.));
}
event measure_tide_gauge (i++) {
  output_gauges (gauges_custom, {eta});
}
/**
~~~gnuplot Tide gauge measurements
reset session
set term svg size 1000,1000
set output "tide_gauges.svg"
set encoding utf8
set autoscale noextend
set tics nomirror
unset key

set multiplot layout 3,1
    unset ylabel
    set margins 8,-1,1,-1
    set title "Left"
    plot "Gauge_left" every 1 using 1:2 with lines title "η"
    set title "Middle"
    plot "Gauge_middle" every 1 using 1:2 with lines title "η"
    set title "Right"
    plot "Gauge_right" every 1 using 1:2 with lines title "η"

    set lmargin -1
    set bmargin -1
    unset title
    set origin 0,-0.005
    set size 1,1
    set border 0
    unset tics
    set ylabel "Free surface elevation (m)" font "," . 14
    set xlabel "Time (s)" font "," . 14
    plot [][0:1] -1     # plot a dummy line out of range
unset multiplot
~~~
*/




/**
Second snapshot.
*/
//event snapshots (t += 10.; t <= 150.) {
event snapshot (t = t1) {
  static int num = 1;
  char filename[128];
  sprintf (filename, "eta%i", num);
  snapshot_fluid (filename, eta, zb, t, ZRANGE_SNAPSHOT);
  num++;
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
event movie (i += animation) {
//event movie (t += 0.3) {
  if (counter == 0) {
    ftime = fopen ("time", "w");
    fanim = fopen ("animation_data", "w");
  }
  foreach(serial)
    fprintf (fanim, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fanim, "\n\n"); 
  fprintf (ftime, "%g\n", t);
  
  counter += 1;
}
#endif







/**
# Final event and output graph of the implementation
We write a file, called `data_out`, to save useful information about the implementation, particularly for the animation script.
*/

event end (t = Tfinal) {
  snapshot_fluid ("eta_END", eta, zb, t, ZRANGE_SNAPSHOT);
  
  #if animation
  fclose (ftime);
  fclose (fanim);
  #endif
  FILE * fp = fopen ("data_out", "w");
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
  fprintf (fp, "counter= %i\n", counter);
  fprintf (fp, "min(zb)= %g\n", statsf(zb).min);
  fclose (fp);

  discrete_Fourier_transform_from_file ("Gauge_left", 0, "Left gauge");
  discrete_Fourier_transform_from_file ("Gauge_middle", 0, "Middle gauge");
  discrete_Fourier_transform_from_file ("Gauge_right", 0, "Right gauge");
}

/** ![Fourier transform of the signal on the left border.](beach_waves/Gauge_left_DFT.pdf) */

/**
Finally we plot the mean and max values of the unknowns ($h$ and $u$)
and information about the time step during the simulation.
*/


/**
~~~gnuplot Values over time
reset session
set term svg size 750,750
#set termoption dashed
set output "log.svg"
set encoding utf8
set key left
C0 = "#9400d3"
C1 = "#009e73"
set multiplot layout 3,1

set xlabel "Time"
set ylabel "Maximum depth" tc rgb C0
set y2label "Mean depth" tc rgb C1
set ytics format "%g" tc rgb C0
set y2tics format "%3.3g" tc rgb C1
plot \
     "log" using 2:4 with lines title "max(h)" lc rgb C0, \
     "log" using 2:5 axes x1y2 with lines title "~h‾ " lc rgb C1

set xlabel "Time"
set ylabel "Maximum velocity" tc rgb C0
set y2label "Mean velocity" tc rgb C1
plot \
     "log" using 2:6 with lines title "max(u)" lc rgb C0, \
     "log" using 2:7 axes x1y2 with lines title "~u‾ " lc rgb C1
     
set xlabel "Iterations"
set ylabel "Time" tc rgb C0
set y2label "Time step" tc rgb C1
set y2tics format "%0.0e" tc rgb C1
set logscale y2 10
plot \
     "log" using 1:2 with lines title "t" lc rgb C0, \
     "log" using 1:3 axes x1y2 with lines title "Δt" lc rgb C1

unset multiplot
~~~
*/
