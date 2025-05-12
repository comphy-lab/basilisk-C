/**
This is a code to produce graphs as post-treatment of the main programm, Port_CFDwavemaker.
*/
#define PATH_TO_RESULTS "../Port_waves/"

/**
We include the same headers as in Port_waves.
*/
#include "grid/multigrid.h"
#include "spherical.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"


#include "terrain.h"
#include "input.h"
#include "output.h"
#include "gauges.h"
#include "view.h"


/**
We include my "custom" headers containing routines to output *via* gnuplot and to discretly Fourier transform.
*/
#include "output_gnuplot.h"
#include "DFT.h"
#include "waves.h"



/**
We define some fields:
a scalar `H` to store the water depth,
a vector `U` to store the average velocities at any given time,
a scalar `Zl` to store the "height" of each layer 
*/

scalar H[];
vector U[];
scalar Zl, volume_l; // Altitude of (the middle of) each layer and volume of each layer;


double box_domain[2][2]; // For my output_gnuplot function
double lat_degree_m, lon_degree_m;

wave_spectrum waves;


double Tfinal;
int counter = 0;
int main ()
{
  Radius = 6371220.; // Earth radius (in meters)

  G = 9.81;
  
  CFL_H = 1.;
  breaking = 0.07;


  /**
     We read the simulation data from the folder with the results from the "original" run.
  */
  char filename[128], buffer[256];
  sprintf (filename, "%sdata_out", PATH_TO_RESULTS);
  fprintf (stderr, "%s\n", filename);
  FILE * fp = fopen (filename, "r");
  fgets (buffer, sizeof(buffer), fp); // Jumping a line
  fscanf (fp, "%lg\n", &X0);
  fprintf (stderr, "X0 = %g\n", X0);
  fgets (buffer, sizeof(buffer), fp);
  fscanf (fp, "%lg\n", &Y0);
  fprintf (stderr, "Y0 = %g\n", Y0);
  fgets (buffer, sizeof(buffer), fp);
  fscanf (fp, "%lg\n", &L0);
  fprintf (stderr, "L0 = %g\n", L0);
  fgets (buffer, sizeof(buffer), fp);
  fscanf (fp, "%i\n", &N);
  fprintf (stderr, "N = %i\n", N);
  fgets (buffer, sizeof(buffer), fp);
  fscanf (fp, "%i\n", &nl);
  fprintf (stderr, "nl = %i\n", nl);
  fgets (buffer, sizeof(buffer), fp);
  fgets (buffer, sizeof(buffer), fp);
  fgets (buffer, sizeof(buffer), fp);
  fscanf (fp, "%lg\n", &Tfinal);
  fprintf (stderr, "Tfinal = %g\n", Tfinal);
  fclose (fp);
  
  box_domain[0][0] = X0; box_domain[0][1] = Y0;
  box_domain[1][0] = X0 + L0; box_domain[1][1] = Y0 + L0;

  lat_degree_m = 110000.; // One degree of latitude in meters
  lon_degree_m = cos((X0+L0)/2.*pi/180.) * 110000.; // One degree of longitude in meters (measured at the center of the domain)

  sprintf (filename , "%sPort_waves.waveinput", PATH_TO_RESULTS);
  waves = init_waves (filename);
  
  run ();
}



/**
Measured fields.
*/
scalar eta_mean[], eta_stddev[], eta_max[]; // Free surface elevation
vector velocity_mean[];
scalar velocity_norm_mean[]; // Velocity
scalar kinetic_energy[], kinetic_energy_mean[]; // Kinetic energy
scalar potential_energy[], potential_energy_mean[]; // Potential energy
scalar total_energy[], total_energy_mean[]; // Total energy




/**
# Initialization
*/

event init (i = 0) {
  srand (time(NULL)); // Random seed
  Zl = new scalar[nl];
  volume_l = new scalar[nl];

  dt = DT;
  /**
     We *resore* the simulation using the dumpfile created at the end of the original run.
   */

  char filename[128];
  sprintf (filename, "%sEND.dump", PATH_TO_RESULTS);
  if (restore (file = filename)) {
    fprintf (stderr, "Restarting simulation from i = %i, going to t = %g\n", i, Tfinal);
  } else {
    fprintf (stderr, "ERROR OCCURED WHEN RESTORING.\n");
  }
  output_ppm (zb, n = 1024, file = "zb.ppm");
}






#define COLOR_BOUNDARIES_ETA min = -0.5, max = 0.5
#define COLOR_BOUNDARIES_ETA_GNU 0., 0.25
/**
# Final event
We plot everything we want to.
*/
event end (t = Tfinal + 0.1) {
//event end (i = 1) {
  scalar mask_dry[];
  foreach()
    mask_dry[] = eta[] - zb[] - 1e-6;
  output_ppm (eta, n = 1024, file = "eta.ppm", mask = mask_dry, COLOR_BOUNDARIES_ETA);
  output_ppm (eta_mean, n = 1024, file = "eta_mean.ppm", mask = mask_dry, COLOR_BOUNDARIES_ETA);
  output_ppm (eta_stddev, n = 1024, file = "eta_stddev.ppm", mask = mask_dry, COLOR_BOUNDARIES_ETA);
  output_ppm (eta_max, n = 1024, file = "eta_max.ppm", mask = mask_dry, COLOR_BOUNDARIES_ETA);
  output_ppm (velocity_mean.x, n = 1024, file = "U.x_mean.ppm", mask = mask_dry);
  output_ppm (velocity_mean.y, n = 1024, file = "U.y_mean.ppm", mask = mask_dry);
  output_ppm (velocity_norm_mean, n = 1024, file = "U_norm_mean.ppm", mask = mask_dry);
  output_ppm (kinetic_energy, n = 1024, file = "ke.ppm", mask = mask_dry);
  output_ppm (kinetic_energy_mean, n = 1024, file = "ke_mean.ppm", mask = mask_dry);
  output_ppm (potential_energy, n = 1024, file = "pe.ppm", mask = mask_dry);
  output_ppm (potential_energy_mean, n = 1024, file = "pe_mean.ppm", mask = mask_dry);
  output_ppm (total_energy, n = 1024, file = "Te.ppm", mask = mask_dry);
  output_ppm (total_energy_mean, n = 1024, file = "Te_mean.ppm", mask = mask_dry);

  output_gnuplot_heatmap (eta, N, "gnuplot_eta", COLOR_BOUNDARIES_ETA_GNU, box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "η (m)", "Free surface elevation");
  output_gnuplot_heatmap (eta_mean, N, "gnuplot_eta_mean", -0.25, 0.25, box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "η (m)", "Mean free surface elevation");
  output_gnuplot_heatmap (eta_stddev, N, "gnuplot_eta_stddev", -0.25, 0.25, box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "η (m)", "Standard deviation of free surface elevation");
  output_gnuplot_heatmap (eta_max, N, "gnuplot_eta_max", COLOR_BOUNDARIES_ETA_GNU, box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "η (m)", "Maximum free surface elevation");
  output_gnuplot_heatmap (velocity_mean.x, N, "gnuplot_U.x_mean", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "U.x (m/s)", "Mean velocity (x)");
  output_gnuplot_heatmap (velocity_mean.y, N, "gnuplot_U.y_mean", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "U.y (m/s)", "Mean velocity (y)");
  output_gnuplot_heatmap (velocity_norm_mean, N, "gnuplot_U_norm_mean", 0., 0., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "||U|| (m/s)", "Mean velocity norm");
  output_gnuplot_heatmap (kinetic_energy, N, "gnuplot_ke", 0., 40., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "ke (J/m^3)", "Kinetic energy");
  output_gnuplot_heatmap (kinetic_energy_mean, N, "gnuplot_ke_mean", 0., 20., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "ke (J/m^3)", "Mean kinetic energy");
  output_gnuplot_heatmap (potential_energy, N, "gnuplot_pe", 0., 50., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "pe (J/m^3)", "Potential energy");
  output_gnuplot_heatmap (potential_energy_mean, N, "gnuplot_pe_mean", 0., 20., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "pe (J/m^3)", "Mean potential energy");
  output_gnuplot_heatmap (total_energy, N, "gnuplot_Te", 0., 70., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "Te (J/m^3)", "Total energy");
  output_gnuplot_heatmap (total_energy_mean, N, "gnuplot_Te_mean", 0., 40., box_domain, 1, mask_dry, "Longitude (°)", "Latitude (°)", "Te (J/m^3)", "Mean total energy");
  
/**
We write some information about the simulation in the file `data_out`.
*/
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

  char filename[128]; // WARNING : the graphs produced by this will be located in the PATH_TO_RESULTS folder
  sprintf (filename, "%sGAUGE_Real_tide_gauge", PATH_TO_RESULTS);
  discrete_Fourier_transform_from_file (filename, 0, "Tide gauge location");
  sprintf (filename, "%sGAUGE_Port_entrance", PATH_TO_RESULTS);
  discrete_Fourier_transform_from_file (filename, 0, "Port entrance");
  sprintf (filename, "%sGAUGE_Bassin_Napoleon", PATH_TO_RESULTS);
  discrete_Fourier_transform_from_file (filename, 0, "Bassin Napoléon");
}


/** ![Maximum wave amplitude.](post_treatment/gnuplot_eta_max.png) */

/** ![Mean total energy.](post_treatment/gnuplot_Te_mean.png) */


