/**
This work is an adaptation of my [first](http://basilisk.fr/sandbox/Gwendal_Leger/wave_topo/beach_topo.c) 
and [second](http://basilisk.fr/sandbox/Gwendal_Leger/waves_1D/waves_1D.c) examples to two dimensions,
based on [shoal-ml.c](http://basilisk.fr/src/examples/shoal-ml.c) and [beach.c](http://basilisk.fr/src/test/beach.c),
using the wave theory based on the relation dispersion derived from Euler's equation, $\omega=\sqrt{gk\tanh(kH)}$.\
The case studied here is the port of Boulogne-sur-mer, located in northern France.
*/

/** ![Animation of the wave elevation on Boulogne-sur-mer. Dark blue is $<-0.5$ meters and Dark red is $>0.5$ meters.](Port_waves/eta.mp4) */

/**
We use the multigrid mesh, with the multi-layer solver.\
We also use the `spherical` module to allow us the study of a fluid at the surface of a sphere (here the Earth)
and the `terrain` module to import custom bathymetry in the form of a kdtree.
*/

#define ML 2


#include "grid/multigrid.h"

#include "spherical.h"
#include "terrain.h"

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"

#include "input.h"
#include "output.h"
#include "gauges.h"
#include "view.h"


/**
We include my "custom" headers containing routines to output *via* gnuplot and to discretly Fourier transform,
as well as a header to produce waves (according to the relation dispersion $\omega=\sqrt{gk\tanh(kH)}$).
*/
#include "output_gnuplot.h"
#include "DFT.h"
#include "waves.h"

/**
We define the level of mesh refinement, $N$ wille be equal to $2^{LEVEL}$ and the final time `Tfinal`.
*/
#define LEVEL 8

#define Tfinal 3600.



/**
Option to save an animation of the simulation,
with the value of `animation` being the number of time steps between frames.\
Multpile outputs are possible, 
`PPM` means using the `output_ppm` option, `anim_DFT` means animating the Fourier transform,
`GNUPLOT` means using gnuplot to plot in 3D
and `GNUPLOT_HEATMAP` produces an animated heatmap.
Note that all these options, except `PPM`, can slow down the code a lot because of the part where it's just gnuplot animating.
*/
#define animation 5
#if animation
#define PPM 1
#define anim_DFT 1 // Can slow down the simulation A LOT
#define GNUPLOT 0 // VERY SLOW
#define GNUPLOT_HEATMAP 0 // SLOW
#if GNUPLOT_HEATMAP
FILE * fhtmp;
#endif
#endif
int counter = 0;


/**
We define some fields:
a scalar `H` to store the water depth,
a vector `U` to store the average velocities at any given time,
a scalar `Zl` to store the "height" of each layer 
and a scalar `volume_l` to store the volume of each cell in each layer.
*/

scalar H[];
vector U[];
scalar Zl, volume_l; // Altitude of (the middle of) each layer and volume of each layer;


double box_domain[2][2]; // For my output_gnuplot function
double lat_degree_m, lon_degree_m;

wave_spectrum waves;


int main ()
{
/**
 We define the domain as a square of length 0.055 degrees ($\approx$ 5 kilometers), centered on the port of Boulogne-sur-mer.
*/
  Radius = 6371220.; // Earth radius (in meters).
  size (0.055); origin (1.6023 - L0, 50.7343 - L0*0.5);

  G = 9.81;
  
  nl = ML;
  CFL_H = 1.;
  breaking = 0.07;
 
  N = 1 << LEVEL;
/**
We force the time step not to be too large in case of lake-at-rest initialization, because the unconditionnal stability of the multilayer solver coupled with a null velocity produces an "inifite" time step.
*/
  DT = 5.;

  box_domain[0][0] = X0; box_domain[0][1] = Y0; // This is for output_gnuplot purposes
  box_domain[1][0] = X0 + L0; box_domain[1][1] = Y0 + L0;

  lat_degree_m = 110000.; // One degree of latitude in meters
  lon_degree_m = cos((X0+L0)/2.*pi/180.) * 110000.; // One degree of longitude in meters (measured at the center of the domain)

  /**
     We setup the wave spectrum using the input data of the file `waves_input.h`.
  */
  waves = init_waves ();
  
  run ();
}


/**
# Initialization
*/
event init (i = 0) {
  Zl = new scalar[nl];
  volume_l = new scalar[nl];

  dt = DT;
  /**
     We use data from the [french Marine's Hydrographic and Oceanic Service (SHOM) database](https://data.shom.fr),
     and we convert it to a kd-tree format using `xyz2kdt` ([documentation on xyz2kdt](https://gfs.sourceforge.net/wiki/index.php/Xyz2kdt)).\
We had to manually complete the seawall because in the data it had many holes due to the data resolution of 10 meters.
   */
  terrain (zb, "Port_waves.BOULOGNE-SUR-MER", "Port_waves.seawall", NULL);
  // Initialization of depth h, velocity u and free surface elevation eta
  foreach() {
    H[] = max (0., -zb[]); // Lake at rest
    double zl = zb[];
    foreach_layer() {
      h[] = H[] / (double)nl;
      zl += h[]*0.5;
      Zl[] = zl;
      zl += h[]*0.5;
    }
    foreach_dimension() {
      foreach_layer() {
	U.x[] = 0.;
	u.x[] = 0.;
      }
    }
    eta[] = zb[] + H[];
  }

  output_ppm (eta, n = 1024, file = "eta0.ppm", min = -1., max = 1.);
  output_ppm (zb, n = 1024, file = "zb.ppm");
    
  /**
     Boundary conditions using the wave theory.
  */
  u.n[left] = VeloX ((x-X0)*lon_degree_m, (y-Y0)*lat_degree_m, Zl[], t, waves, H[]);
  u.t[left] = VeloY ((x-X0)*lon_degree_m, (y-Y0)*lat_degree_m, Zl[], t, waves, H[]);
  w[left] = VeloZ ((x-X0)*lon_degree_m, (y-Y0)*lat_degree_m, Zl[]+h[]*0.5, t, waves, H[]);
  //u.n[left] = -radiation (0.);
  u.n[right] = radiation (0.);
  u.n[top] = radiation (0.);
  u.n[bottom] = -radiation (0.);
}




/**
# Events
*/
/**
This is where we compute the water depth `H` and the (vertical) mean velocity `U`.
*/
event compute_vertical_scalars (i++) {
  foreach() {
    H[] = eta[] - zb[];
    foreach_dimension() {
      U.x[] = 0.;
      if (H[] > dry) {
	foreach_layer()
	  U.x[] += u.x[]*h[];
	U.x[] /= H[];
      }
    }
    double zl = zb[];
    foreach_layer() {
      zl += 0.5*h[];
      Zl[] = zl;
      volume_l[] = dv()*h[];
      zl += 0.5*h[];
    }
  }
}


/**
This is the friction event from [beach.c](http://basilisk.fr/src/test/beach.c), adapted to a 2D case.
*/
#if 1
event friction (i++) {
  foreach() {
    double Qx = H[]*U.x[];
    double Qy = H[]*U.y[];
    double a = (H[] < dry ? HUGE : 1. + 5e-3*dt*sqrt(sq(Qx) + sq(Qy))/sq(H[]));
    foreach_dimension() {
      foreach_layer()
	u.x[] /= a;
      U.x[] = 0.;
      if (H[] > dry) { // If the cell isn't dry, we need to re-compute U
	foreach_layer()
	  U.x[] += u.x[]*h[];
	U.x[] /= H[];
      }
    }
  }
}
#endif // if 1

/**
We measure different fields, the water depth average, standard deviation and maximum value in time.
*/
scalar eta_mean[], eta_stddev[], eta_max[]; // Free surface elevation
vector velocity_mean[];
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
    foreach_dimension()
      velocity_mean.x[] += U.x[]*dt;
    velocity_norm_mean[] = sqrt (sq(U.x[]) + sq(U.y[]))*dt;
    kinetic_energy[] = 0.;
    foreach_layer () 
      kinetic_energy[] += 0.5*volume_l[]*(sq(u.x[]) + sq(u.y[]) + sq(w[]));
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
  stats stats_ux = statsf (U.x);
  stats stats_uy = statsf (U.y);
  double mean_ux = stats_ux.sum/stats_ux.volume;
  double mean_uy = stats_uy.sum/stats_uy.volume;
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
    foreach_layer()
      sum_volume += volume_l[];
  }
  if (i == 0) {
    fprintf (stderr, "log  i     t   dt        h: mean; max          eta: min; mean; max     u: mean; max; mean; max     ke: mean; max    pe: mean; max    Te: mean; max    E:  k   p   T    volume/volume0\n");
    volume_t0 = sum_volume;
  }
  fprintf (stderr, "log %i %5.5g %3.3e   %g %g   %g %g %g  %g %g %g %g   %g %g   %g %g   %g %g   %g %g %g   %g\n",
	   i, t, dt,
	   mean_h, stats_h.max,
	   stats_eta.min, mean_eta, stats_eta.max,
	   mean_ux, stats_ux.max, mean_uy, stats_uy.max,
	   mean_ke, stats_ke.max,
	   mean_pe, stats_pe.max,
	   mean_Te, stats_Te.max,
	   Ek, Ep, ET,
	   (sum_volume - volume_t0)/volume_t0);
}



/**
We use the `gauges.h` module to measure the evolution of the surface in different places, 
at the location of a real tide gauge in Boulogne-sur-mer (to compare later),
at the entrance of the port (le Caisson)
and at the back of the port, near the bassin Napoléon.
*/
const int n_gauges = 3; // Number of tide gauges
Gauge gauge_boulogne[] = {
			  {"GAUGE_Real_tide_gauge", 1.57688, 50.72738, "Boulogne-sur-mer tide gauge (real location)"}, // The actual tide gauge is on the shore so we shifted its position by a 1/1000th of degree west.
			  {"GAUGE_Port_entrance", 1.5699, 50.7453, "Port entrance tide gauge (Le Caisson)"},
			  {"GAUGE_Bassin_Napoleon", 1.6007, 50.7236, "Back of bassin Napoléon tide gauge"},
			  {NULL}};
event init_gauges (i = 0) {
  FILE * fp = fopen ("tide_gauges_positions", "w");
  for (int ii = 0; ii < n_gauges; ii++)
    fprintf (fp, "%g %g\n", gauge_boulogne[ii].x, gauge_boulogne[ii].y);
  fclose (fp);
  output_gnuplot_heatmap (zb, N, "tide_gauges_locations", 0., 0., box_domain, 0, zb, "Longitude (°)", "Latitude (°)", "zb (m)", "Tide gauges locations");
  fp = fopen ("plot.plot", "a");
  fprintf (fp, ", 'tide_gauges_positions' using 1:2:0 with points lc rgb '#FFFFFF'\n");
  fclose (fp);
  system ("gnuplot plot.plot");
}
event measure_tide_gauge (i++) {
  output_gauges (gauge_boulogne, {eta});
}
/**
~~~gnuplot Tide gauge measurements
reset session
set term svg size 750,750
set output "tide_gauges_b-s-m.svg"
set encoding utf8
set autoscale noextend
unset key

set multiplot layout 3,1
    unset ylabel
    set margins 8,-1,1,-1
    set title "Real tide gauge"
    plot "GAUGE_Real_tide_gauge" every 1 using 1:2 with lines title "η"
    set title "Port entrance (Le Caisson)"
    plot "GAUGE_Port_entrance" every 1 using 1:2 with lines title "η"
    set title "Bassin Napoléon"
    plot "GAUGE_Bassin_Napoleon" every 1 using 1:2 with lines title "η"

    set lmargin -1
    set bmargin -1
    unset title
    set origin 0,0
    set size 1,1
    set border 0
    unset tics
    set ylabel "Free surface elevation (m)" font "," . 14
    set xlabel "Time (s)" font "," . 14
    plot [][0:1] -1     # plot a dummy line out of range
unset multiplot
~~~
*/




#define COLOR_BOUNDARIES_ETA min = -0.5, max = 0.5
#if animation
FILE * ftime;
event movie_init (i = 0) {
  ftime = fopen ("time", "w");
  for (int ii = 0; ii < n_gauges; ii++) { // We clear the files used for the DFT animation, in case there are still present
    char filename[128];
    sprintf (filename, "%s_animation_data", gauge_boulogne[ii].name);
    FILE * fp = fopen (filename, "w");
    fclose (fp);
    sprintf (filename, "%s_animation_time", gauge_boulogne[ii].name);
    fp = fopen (filename, "w");
    fclose (fp);
  }
}
event movie (i += animation) {
//event movie (t += 0.1) {
  counter += 1;
  
#if PPM
  scalar mask_dry[];
  foreach()
    mask_dry[] = eta[] - zb[] - 1e-6;
  output_ppm (eta, n = 1024, file = "eta.mp4", mask = mask_dry, COLOR_BOUNDARIES_ETA);
  output_ppm (eta, n = 1024, file = "eta_zoom.mp4", box = {{1.5623, 50.7165}, {1.602, 50.732}}, min = -0.1, max = 0.1, mask = mask_dry);
  output_ppm (U.x, n = 1024, file = "U_x.mp4", mask = mask_dry);
  output_ppm (U.y, n = 1024, file = "U_y.mp4", mask = mask_dry);
  output_ppm (w, n = 1024, file = "w.mp4", mask = mask_dry);
  output_ppm (kinetic_energy, n = 1024, file = "ke.mp4", mask = mask_dry);
  output_ppm (potential_energy, n = 1024, file = "pe.mp4", mask = mask_dry);
  output_ppm (total_energy, n = 1024, file = "Te.mp4", mask = mask_dry);
  scalar velocity_norm[];
  foreach()
    velocity_norm[] = sqrt (sq(U.x[]) + sq(U.y[]));
  output_ppm (velocity_norm, n = 1024, file = "U_norm.mp4", mask = mask_dry);
  scalar actual_eta_mean[], actual_eta_stddev[];
  vector actual_velocity_mean[];
  scalar actual_velocity_norm_mean[];
  scalar actual_kinetic_energy_mean[];
  scalar actual_potential_energy_mean[];
  scalar actual_total_energy_mean[];
  if (t > 0.) {
    const double inv_t = 1./t;
    foreach() {
      actual_eta_mean[] = eta_mean[]*inv_t;
      actual_eta_stddev[] = sqrt (eta_stddev[]*inv_t);
      foreach_dimension()
	actual_velocity_mean.x[] = velocity_mean.x[]*inv_t;
      actual_velocity_norm_mean[] = velocity_norm_mean[]*inv_t;
      actual_kinetic_energy_mean[] = kinetic_energy_mean[]*inv_t;
      actual_potential_energy_mean[] = potential_energy_mean[]*inv_t;
      actual_total_energy_mean[] = total_energy_mean[]*inv_t;
    }
    output_ppm (actual_eta_mean, n = 1024, file = "cvgce_eta_mean.mp4", mask = mask_dry, COLOR_BOUNDARIES_ETA);
    output_ppm (actual_eta_stddev, n = 1024, file = "cvgce_eta_stddev.mp4", mask = mask_dry, COLOR_BOUNDARIES_ETA);
    output_ppm (eta_max, n = 1024, file = "cvgce_eta_max.mp4", mask = mask_dry, COLOR_BOUNDARIES_ETA);
    output_ppm (actual_velocity_mean.x, n = 1024, file = "cvgce_U.x_mean.mp4", mask = mask_dry);
    output_ppm (actual_velocity_mean.y, n = 1024, file = "cvgce_U.y_mean.mp4", mask = mask_dry);
    output_ppm (actual_velocity_norm_mean, n = 1024, file = "cvgce_U_norm_mean.mp4", mask = mask_dry);
    output_ppm (actual_kinetic_energy_mean, n = 1024, file = "cvgce_ke_mean.mp4", mask = mask_dry);
    output_ppm (actual_potential_energy_mean, n = 1024, file = "cvgce_pe_mean.mp4", mask = mask_dry);
    output_ppm (actual_total_energy_mean, n = 1024, file = "cvgce_Te_mean.mp4", mask = mask_dry);
  }
  #endif

  #if GNUPLOT || GNUPLOT_HEATMAP
  char filename[32];
  #endif
  
  #if GNUPLOT
  sprintf (filename, "testtopo_%5.5i.txt", counter);
  FILE * fp = fopen (filename, "w");
  foreach(serial)
    fprintf (fp, "%g %g %g\n", x, y, eta[]);
  fclose (fp);
  #endif

  #if GNUPLOT_HEATMAP
  int ii = 0;
  if (counter == 1)
    fhtmp = fopen ("topo2D_heatmap_data", "w");
  foreach(serial) {
    ii += 1;
      fprintf (fhtmp, "%g %g %g\n", x, y, eta[]);
    if (ii % N == 0)
      fprintf (fhtmp, "\n");
  }
  fprintf (fhtmp, "\n\n");
  #endif
  
  fprintf (ftime, "%g\n", t);

  #if anim_DFT
  discrete_Fourier_transform_animated (gauge_boulogne[0].name, i + 1, t);
  discrete_Fourier_transform_animated (gauge_boulogne[1].name, i + 1, t);
  discrete_Fourier_transform_animated (gauge_boulogne[2].name, i + 1, t);
  #endif
}
#endif



/**
# Final event
*/

event end (t = Tfinal) {
//event end (i = 1) {
  foreach() {
    eta_mean[] /= t;
    eta_stddev[] = sqrt (eta_stddev[]/t);
    foreach_dimension()
      velocity_mean.x[] /= t;
    velocity_norm_mean[] /= t;
    kinetic_energy_mean[] /= t;
    potential_energy_mean[] /= t;
    total_energy_mean[] /= t;
  }
  scalar mask_dry[];
  foreach()
    mask_dry[] = eta[] - zb[] - 1e-6;
  output_ppm (eta, n = 1024, file= "eta.ppm", mask = mask_dry, COLOR_BOUNDARIES_ETA);
  output_ppm (eta_mean, n = 1024, file= "eta_mean.ppm", mask = mask_dry, COLOR_BOUNDARIES_ETA);
  output_ppm (eta_stddev, n = 1024, file= "eta_stddev.ppm", mask = mask_dry, COLOR_BOUNDARIES_ETA);
  output_ppm (eta_max, n = 1024, file= "eta_max.ppm", mask = mask_dry, COLOR_BOUNDARIES_ETA);
  output_ppm (velocity_mean.x, n = 1024, file= "U.x_mean.ppm", mask = mask_dry);
  output_ppm (velocity_mean.y, n = 1024, file= "U.y_mean.ppm", mask = mask_dry);
  output_ppm (velocity_norm_mean, n = 1024, file= "U_norm_mean.ppm", mask = mask_dry);
  output_ppm (kinetic_energy, n = 1024, file = "ke.ppm", mask = mask_dry);
  output_ppm (kinetic_energy_mean, n = 1024, file = "ke_mean.ppm", mask = mask_dry);
  output_ppm (potential_energy, n = 1024, file = "pe.ppm", mask = mask_dry);
  output_ppm (potential_energy_mean, n = 1024, file = "pe_mean.ppm", mask = mask_dry);
  output_ppm (total_energy, n = 1024, file = "Te.ppm", mask = mask_dry);
  output_ppm (total_energy_mean, n = 1024, file = "Te_mean.ppm", mask = mask_dry);

  output_gnuplot_heatmap (eta, N, "gnuplot_eta", -0.5, 0.5, box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "eta", "eta");
  output_gnuplot_heatmap (eta_mean, N, "gnuplot_eta_mean", -0.5, 0.5, box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "eta", "eta (mean)");
  output_gnuplot_heatmap (eta_stddev, N, "gnuplot_eta_stddev", -0.5, 0.5, box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "eta", "eta standard deviation");
  output_gnuplot_heatmap (eta_max, N, "gnuplot_eta_max", -0.5, 0.5, box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "eta", "eta maximum");
  output_gnuplot_heatmap (velocity_mean.x, N, "gnuplot_U.x_mean", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "U.x", "Velocity (x) mean");
  output_gnuplot_heatmap (velocity_mean.y, N, "gnuplot_U.y_mean", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "U.y", "Velocity (y) mean");
  output_gnuplot_heatmap (kinetic_energy, N, "gnuplot_ke", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "ke", "Kinetic energy");
  output_gnuplot_heatmap (kinetic_energy_mean, N, "gnuplot_ke_mean", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "ke", "Kinetic energy mean");
  output_gnuplot_heatmap (potential_energy, N, "gnuplot_pe", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "pe", "Potential energy");
  output_gnuplot_heatmap (potential_energy_mean, N, "gnuplot_pe_mean", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "pe", "Potential energy");
  output_gnuplot_heatmap (total_energy, N, "gnuplot_Te", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "Te", "Total energy");
  output_gnuplot_heatmap (total_energy_mean, N, "gnuplot_Te_mean", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "Te", "Total energy mean");
  
/**
We write some information about the simulation in the file `data_out`, which is particularly useful for the animation script.
*/
  dump ("END.dump");
  FILE * fp = fopen ("data_out", "w");
  // "Physical parameters"
  fprintf (fp, "X0\n%g\n", X0);
  fprintf (fp, "Y0\n%g\n", Y0);
  fprintf (fp, "L0\n%g\n", L0);
  // Numerical parameters
  fprintf (fp, "N\n%i\n", N);
  fprintf (fp, "nl\n%i\n", nl);
  // End state
  fprintf (fp, "i\n%i\n", i);
  fprintf (fp, "t\n%g\n", t);
  // Animation data
  fprintf (fp, "counter\n%i\n", counter);
  fprintf (fp, "min(zb)\n%g\n", statsf(zb).min);
  fclose (fp);

  
  #if animation
  fclose (ftime);
  #if GNUPLOT || GNUPLOT_HEATMAP
  FILE * fplot;
  #endif
  #if GNUPLOT
  fplot = fopen ("plot_anim_topo2D.plot", "w");
  fprintf (fplot,
  	   "set term gif animate delay 5\n"
  	   "set output 'topo2D.gif'\n"
  	   "unset key\n"
  	   "set xr[%g:%g]\n"
  	   "set yr[%g:%g]\n"
  	   "set zr[%g:%g]\n"
  	   "ftime = 'time'\n"
  	   "getValue (row, col, filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')\n"
  	   "do for [i=1:%d] {\n"
  	   "\tset title getValue(i,1,ftime)\n"
  	   "\tfile = sprintf('testtopo_%%05.0f.txt', i)\n"
  	   "#\tset view 60, i%%360, 1, 1\n"
  	   "\tsplot 'bathymetry_position' using 1:2:3 with points pointtype 0 lc rgb '#000000',"
  	   "file using 1:2:3 with points pointtype 0 lc rgb '#1234FF'\n"
  	   "}\n",
  	   X0, X0+L0, Y0, Y0+L0, statsf(zb).min, statsf(eta_max).max, counter);
  fclose (fplot);
  system ("gnuplot plot_anim_topo2D.plot");
  #endif
  #if GNUPLOT_HEATMAP
  fplot = fopen ("plot_anim_topo2D_heatmap.plot", "w");
  fprintf (fplot,
	   "set term gif animate delay 5\n"
	   "set output 'topo2D_heatmap.gif'\n"
	   "set view map scale 1\n"
	   "set size square\n"
	   "set pm3d implicit\n"
	   "unset key\n"
	   "set xr[%g:%g]\n"
	   "set yr[%g:%g]\n"
	   "set zr[%g:%g]\n"
	   "# set palette defined (0 '#32369c', 1 '#0098fe', 2 '#05cd67', 3 '#5ddf79', 4 '#c9f48e', 5 '#e4dc8a', 6 '#aa926b')\n"
	   "ftime = 'time'\n"
	   "file = 'topo2D_heatmap_data'\n"
	   "getValue (row, col, filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')\n"
	   "do for [i=1:%d] {\n"
	   "\tset title getValue(i,1,ftime)\n"
	   "\tsplot file index (i-1) using 1:2:3 with pm3d\n"
	   "}\n",
	   X0, X0+L0, Y0, Y0+L0, statsf(zb).min, statsf(eta_max).max, counter);
  fclose (fplot);
  system ("gnuplot plot_anim_topo2D_heatmap.plot");
  #endif
  #endif

  discrete_Fourier_transform_from_file (gauge_boulogne[0].name, 0, "Real tide gauge");
  discrete_Fourier_transform_from_file (gauge_boulogne[1].name, 0, "Port entrance");
  discrete_Fourier_transform_from_file (gauge_boulogne[2].name, 0, "Bassin Napoléon");
  #if anim_DFT
  discrete_Fourier_transform_animated (gauge_boulogne[0].name, 0, t);
  discrete_Fourier_transform_animated (gauge_boulogne[1].name, 0, t);
  discrete_Fourier_transform_animated (gauge_boulogne[2].name, 0, t);
  #endif  
}


/**
# Some results
*/

/** ![Maximum wave amplitude.](Port_waves/eta_max.ppm) */

/** ![Animation of the Fourier transform of the signal at the back of the port, near bassin Napoléon.
    Note the fundamental frequence being excited.](Port_waves/GAUGE_Bassin_Napoleon_DFT.gif) */

/** ![Fourier transform of the signal at the entrance of the port.
    Note that we see much more frequencies than at the back of the port.](Port_waves/GAUGE_Port_entrance_DFT.gif) */



 



/**
~~~gnuplot Values over time (note that there are two scales !).
reset session
set term svg size 1200,1200
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
     "log" using 3:11 with lines title "M(u_x)" dashtype 3 lc rgb C0, \
     "log" using 3:10 axes x1y2 with lines title "~u‾ _x " lc rgb C1, \
     "log" using 3:13 with lines title "M(u_y)" dashtype 2 lc rgb C0, \
     "log" using 3:12 axes x1y2 with lines title "~u‾ _y " dashtype 4 lc rgb C1

unset logscale x
set xlabel "Time (s)"
set ylabel "Energy (J)" tc rgb C0
set ytics format "%0.0e" nomirror
unset logscale y
set y2label "Relative volume variation" tc rgb C1
set y2tics format "%2.2g" tc rgb C1
plot \
     "log" using 3:14 with lines title "ke" lc rgb C0, \
     "log" using 3:16 with lines title "pe" dashtype 3 lc rgb C0, \
     "log" using 3:($18*0.5) with lines title "Te/2" dashtype 4 lc rgb C0, \
     "log" using 3:23 axes x1y2 with lines title "(V-V0)/V0" lc rgb C1
     
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
