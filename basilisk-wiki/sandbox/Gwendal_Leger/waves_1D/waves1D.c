/**
This work is an adaptation of my [first example](http://basilisk.fr/sandbox/Gwendal_Leger/wave_topo/beach_topo.c) based on [beach.c](http://basilisk.fr/src/test/beach.c), 
but with wave theory based on the dispersion relation $\omega=\sqrt{gk\tanh(kH)}$.

See [the results for the same code but with the Saint-Venant solver](http://basilisk.fr/sandbox/Gwendal_Leger/waves_1D/waves_1D_SV.c) to see why dispersivity (wave velocity depending on wavelength) is important !
*/

/**
~~~gnuplot Animation
reset session
set term gif animate delay 5
set output 'animation.gif'
unset key
getValue (row, col, filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')
data_out = "data_out"
set xr[getValue (1, 2, data_out):getValue (1, 2, data_out)+getValue (2, 2, data_out)] # X0 and L0 are, respectively, at rows 1 and 2, columns 2 of "data_out"
#set yr[getValue (8, 2, data_out):getValue (9, 2, data_out)]
set yrange [-1.:1.]
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
We define the number of layer, and if that number is 0, we use the Saint-Venant solver instead of the (Euler) multilayer solver.
*/
#define ML 2

#include "grid/multigrid1D.h"

#if ML
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#else
#include "saint-venant.h"
#endif

#include "output.h"
#include "gauges.h"


/**
We include my "custom" headers containing routines to output *via* gnuplot, to discretly Fourier transform and to use wave theory.\
To input wave data, we need to write it directly in the `waves_input.h` file.
*/
#include "output_gnuplot1D.h"
#include "DFT.h"
#include "waves.h"

/**
We define the level of mesh refinement, $N$ wille be equal to $2^{LEVEL}$ and the final time `Tfinal`.
*/
#define LEVEL 10

#define Tfinal 50.



/**
Options to display an animation of the simulation while running and/or to save an animation of the simulation,
with the value of `animation` being the number of time steps between frames.
*/
#define animation_while_running 0

#define animation 5
int counter = 0;
#if animation
FILE *ftime, *fanim;
double maxeta_global = 0.;
#endif


/**
We define some fields:
a scalar `H` to store the water depth,
a vector `U` to store the average velocities at any given time,
a scalar `Zl` to store the "height" (middle) of each layer 
*/

scalar H[], U[];
#if ML
scalar Zl, volume_l; // Altitude of (the middle of) each layer and volume of each layer
#else
scalar Zl[], volume_l[];
#endif

wave_spectrum waves;


int main ()
{
/**
 We define the domain as an interval of length 100 centered on 0.
*/
  size (100.);
  origin (-L0*0.5);
  
  G = 9.81;

  #if ML
  nl = ML;
  CFL_H = 1.;
  breaking = 0.07;
  #else
  nl = 1;
  #endif
 
  N = 1 << LEVEL;
/**
We force the time step not to be too large in case of lake-at-rest initialization, because the unconditionnal stability of the multilayer solver coupled with a null velocity produces an "inifite" time step.
*/
  DT = 0.1;
/**
We setup the wave spectrum that will be used, by using the data from `waves_input.h` file.
*/     
  waves = init_waves ();
  
  run ();
}

/**
# Initialization
*/

event init (i = 0) {
  #if ML
  Zl = new scalar[nl];
  volume_l = new scalar[nl];
  #endif
  dt = DT;
  
  // Initialization of depth h, velocity u and free surface elevation eta
  foreach() {
    zb[] = -10.;
    eta[] = 0.;
    H[] = max (0., eta[] - zb[]); // Lake at rest
    U[] = 0.;
    double zl = zb[];
    #if ML
    foreach_layer()
    #endif
      {
	h[] = H[]/(double)nl;
	zl += h[]*0.5;
	Zl[] = zl;
	u.x[] = 0.;
	#if ML
	w[] = 0.;
	#endif
	zl += h[]*0.5;
      }
  }      
  /**
     Boundary conditions.						\
     We have to be careful because $u_k$ is defined at the middle of layer $k$ while $w$ is defined at the top of the layer $k$.
  */
  u.n[left] = VeloX (x, 0., Zl[] - zb[], t, waves, -zb[]);
  #if ML // We differentiate the boundary conditions so that the fluid "gets out" of the domain easily
  u.n[right] = radiation (0.);
  #else
  u.n[right] = neumann (0.);
  #endif

  snapshot_fluid ("eta0", eta, zb, 0., -1., 1.);
}



/**
# Events
*/
/**
This is where we compute the water depth `H` and the (vertical) mean velocity `U`, as well as "tool" scalars containaing the volume (`volume_l`) and the average altitude (`Zl`) of each layer.
*/
event compute_vertical_scalars (i++) {
  foreach() {
    H[] = eta[] - zb[];
    U[] = 0.;
    if (H[] > dry) {
      #if ML
      foreach_layer()
      #endif
	U[] += u.x[]*h[];
      U[] /= H[];
    }
    double zl = zb[];
    #if ML
    foreach_layer()
    #endif
      {
	zl += 0.5*h[];
	Zl[] = zl;
	zl += 0.5*h[];
	volume_l[] = dv()*h[];
      }
  }
}


/**
This is the friction event from [beach.c](http://basilisk.fr/src/test/beach.c).
*/
#if 1
event friction (i++) {
  foreach() {
    #if ML
    double Q = H[]*U[];
    double a = (H[] < dry ? HUGE : 1. + 5e-3*dt*fabs(Q)/sq(H[]));
    foreach_layer()
      u.x[] /= a;
    U[] = 0.;
    if (H[] > dry) { // If the cell isn't dry, we need to re-compute U
      foreach_layer()
	U[] += u.x[]*h[];
      U[] /= H[];
    }
    #else
    double a = h[] < dry ? HUGE : 1. + 5e-3*dt*norm(u)/h[];
    u.x[] /= a;
    U[] = u.x[];
    #endif
  }
}
#endif // if 1


/**
We measure different fields, the water depth average, standard deviation and maximum value in time.
*/
scalar eta_mean[], eta_stddev[], eta_max[]; // Free surface elevation
scalar velocity_mean[];
scalar velocity_norm_mean[]; // Velocity
scalar kinetic_energy[], kinetic_energy_mean[]; // Kinetic energy
scalar potential_energy[], potential_energy_mean[]; // Potential energy
scalar total_energy[], total_energy_mean[]; // Total energy

event measure_fields (i++) {
  foreach() {
    eta_mean[] += eta[]*dt;
    eta_stddev[] += sq(eta[]*dt - eta_mean[]/(t > 0. ? t : HUGE)); // Not to divide by 0
    if (eta[] > eta_max[])
      eta_max[] = eta[];
    velocity_mean[] += U[]*dt;
    velocity_norm_mean[] = fabs(U[])*dt;
    kinetic_energy[] = 0.;
    #if ML
    foreach_layer ()
    #endif
      kinetic_energy[] += 0.5*volume_l[]*sq(u.x[]);
    potential_energy[] = H[] > dry ? G*sq(eta[])*0.5*dv() : 0.;
    total_energy[] = kinetic_energy[] + potential_energy[];
    kinetic_energy_mean[] += kinetic_energy[]*dt;
    potential_energy_mean[] += potential_energy[]*dt;
    total_energy_mean[] += total_energy[]*dt;
  }
}


double volume_t0;
event logfile (i += 10) {
  //const double one_over_NN = 1./(double)(N*N);
  stats stats_h = statsf (H);
  double mean_h = stats_h.sum/stats_h.volume;
  stats stats_eta = statsf (eta);
  double mean_eta = stats_eta.sum/stats_eta.volume;
  stats stats_u = statsf (U);
  double mean_u = stats_u.sum/stats_u.volume;
  stats stats_ke = statsf (kinetic_energy);
  double mean_ke = stats_ke.sum/stats_ke.volume;
  stats stats_pe = statsf (potential_energy);
  double mean_pe = stats_pe.sum/stats_pe.volume;
  stats stats_Te = statsf (total_energy);
  double mean_Te = stats_Te.sum/stats_Te.volume;
  double Ek = 0., Ep = 0., ET = 0.;
  double sum_volume = 0.;
  foreach (reduction(+:Ek) reduction(+:Ep) reduction(+:ET) reduction(+:sum_volume)) {
    Ek += kinetic_energy[];
    Ep += potential_energy[];
    ET += (kinetic_energy[] + potential_energy[])*0.5;
    #if ML
    foreach_layer()
    #endif
      sum_volume += volume_l[];
  }
  if (i == 0) {
    fprintf (stderr, "log  i     t   dt        h: mean; max          eta: min; mean; max     u: mean; max     ke: mean; max    pe: mean; max    Te: mean; max    E:  k   p   T    volume/volume0\n");
    volume_t0 = sum_volume;
  }
  fprintf (stderr, "log %i %5.5g %3.3e   %g %g   %g %g %g  %g %g   %g %g   %g %g   %g %g   %g %g %g   %g\n",
	   i, t, dt,
	   mean_h, stats_h.max,
	   stats_eta.min, mean_eta, stats_eta.max,
	   mean_u, stats_u.max,
	   mean_ke, stats_ke.max,
	   mean_pe, stats_pe.max,
	   mean_Te, stats_Te.max,
	   Ek, Ep, ET,
	   (sum_volume - volume_t0)/volume_t0);
}


/**
Event from `beach.c` to display an animation of the wave while the programm is running.
*/
#if animation_while_running
event running_movie (i += animation_while_running) {
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



#if animation
event movie (i += animation) {
//event movie (t += 0.1) {
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
# Final event
*/
event end (t = Tfinal) {
//event end (i = 1) {
  foreach() {
    eta_mean[] /= t;
    eta_stddev[] = sqrt (eta_stddev[]/t);
    velocity_mean[] /= t;
    velocity_norm_mean[] /= t;
    kinetic_energy_mean[] /= t;
    potential_energy_mean[] /= t;
    total_energy_mean[] /= t;
  }
  
/**
We write some information about the simulation in the file `data_out`, which is particularly useful for the animation script, and do the Fourier transforms of the signals measured by the "tide gauges".
*/
  dump ("END.dump");
  FILE * fp = fopen ("data_out", "w");
  // "Physical parameters"
  fprintf (fp, "X0= %g\n", X0);
  fprintf (fp, "L0= %g\n", L0);
  // Numerical parameters
  fprintf (fp, "N= %i\n", N);
  fprintf (fp, "nl= %i\n", nl);
  // End state
  fprintf (fp, "i= %i\n", i);
  fprintf (fp, "t= %g\n", t);
  // Animation data
  fprintf (fp, "counter= %i\n", counter);
  fprintf (fp, "min(zb)= %g\n", statsf(zb).min);
  fprintf (fp, "max(eta)= %g\n", statsf(eta_max).max);
  fclose (fp);

  
  #if animation
  fclose (ftime);
  fclose (fanim);
  #endif

  discrete_Fourier_transform_from_file ("Gauge_left", 0, "Left gauge");
  discrete_Fourier_transform_from_file ("Gauge_middle", 0, "Middle gauge");
  discrete_Fourier_transform_from_file ("Gauge_right", 0, "Right gauge");
}





 



/**
~~~gnuplot Values over time
reset session
set term svg size 1000,1000 #size 3820,2150
set termoption dashed
set output "log.svg"
set encoding utf8
set key left
C0 = "#9400d3"
C1 = "#009e73"
C2 = "#ff0000"
set autoscale noextend
set multiplot layout 4,1

set xlabel "Time (s)"
set ylabel "Maximum depth (m)" tc rgb C0
set y2label "Mean depth (m)" tc rgb C1
set ytics format "%g" tc rgb C0 nomirror
set y2tics format "%g" tc rgb C1
plot \
     "log" using 3:6 with lines title "max(h)" lc rgb C0, \
     "log" using 3:5 axes x1y2 with lines title "~h‾ " lc rgb C1

set xlabel "Time (s)"
set ylabel "Maximum velocity (m.s^{-1})" tc rgb C0
set ytics nomirror
set y2label "Mean velocity (m.s^{-1})" tc rgb C1
plot \
     "log" using 3:11 with lines title "M(u_x)" lc rgb C0, \
     "log" using 3:10 axes x1y2 with lines title "~u‾ _x " lc rgb C1

unset logscale x
set xlabel "Time (s)"
set ylabel "Energy (J)" tc rgb C0
set ytics format "%0.0e" nomirror
unset logscale y
set y2label "Relative volume variation" tc rgb C1
set y2tics format "%2.2g" tc rgb C1
plot \
     "log" using 3:12 with lines title "ke" lc rgb C0, \
     "log" using 3:14 with lines title "pe" dashtype 3 lc rgb C0, \
     "log" using 3:($16*0.5) with lines title "Te/2" dashtype 4 lc rgb C0, \
     "log" using 3:21 axes x1y2 with lines title "(V-V0)/V0" lc rgb C1
     
set xlabel "Iterations"
unset logscale y
set ylabel "Time (s)" tc rgb C0
set ytics format "%g" nomirror
set y2label "Time step (s)" tc rgb C1
set y2tics format "%0.0e" tc rgb C1
set logscale y2 10
plot \
     "log" using 2:3 with lines title "t" lc rgb C0, \
     "log" using 2:4 axes x1y2 with lines title "Δt" lc rgb C1

unset multiplot
~~~
*/
