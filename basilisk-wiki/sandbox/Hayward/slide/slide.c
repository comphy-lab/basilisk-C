/**
# Tsunami generation and run-up by a submerging variable density flow

This is lightly adpated from [Lily's version](/sandbox/lbattershill/multilayer-slide/slide.c) to refine the geometry, and add a run-up slope.

*/


#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#define drho(T) T
#include "layered/dr.h"
#include "layered/remap.h"
#include "layered/rpe.h"
#include "layered/perfs.h"

#define nlayers 10                        // Number of layers in simulation
#define resolution 16.                    // Horzontal resoltion (cells per x)

// Densities
double rho_0 = 0.;                        // Water
double rho_1 = 0.01;                      // 
//double rho_2 = -0.01;                     // 





/**
 * Define the geometry. 
*/


#define waterdepth 1.0                    // Bottom water depth
#define slope_exposed_height 2.0          // Height of inital ramp at the foot of the granular flow
#define theta0 0.30                       // Angle of the initial ramp
#define width 1.28                        // hoirzontal width of granular fluid
#define height_cm 0.90                    // Vertical height of granular fluid
#define length_until_runup 25.            // x-position of opposite 'shoreline'
#define theta1 0.25                       // Angle of run-up slope
#define runup_slope_height 3*waterdepth   // Height of run-up slope (to determin domain extent)



#define gradient_m -1.*tan(theta0)                                          // (Calculated) Initial ramp gradient
#define c_intercept (tan(theta0)*width)+slope_exposed_height                // (Calculated) height of initial ramp at x=0
#define domain_length length_until_runup + (runup_slope_height/tan(theta1)) // (Calculated) length of domain





/**
 * Simulation setup. 
*/

double max_T = 0.01;
double min_T = -0.01;
double max_runup_dist = 0.;
double runup_tolerance = 0.001;  // Minimum water height to register as inundated on run-up slope.

int main()
{
  size (domain_length);
  breaking = 0.7;
  N = resolution*domain_length;
  nl = nlayers;
  nu = 0.0001;   // Vertical viscous term
  DT = 0.0025;   // Set time step
  G = 9.81;      // Gravitational acceleration
  cell_lim = mono_limit;
  dry = 0.00002;

  // Make directories
  system ("rm -r gnuplot");
  system ("mkdir gnuplot");
  system ("mkdir gnuplot/all");
  system ("mkdir gnuplot/all/ux");
  system ("mkdir gnuplot/closeup");

  // Free-slip on bottom layer boundary
  const vector slip[] = {HUGE};
  lambda_b = slip;
  
  // Start main loop
  run();
}







/**
 * Initialise domain.
*/

event init (i = 0)
{
  foreach(reduction(max:max_T) reduction(min:min_T)) {
    // BATHMETRY
    zb[] = x < width + (waterdepth+slope_exposed_height)/tan(theta0) ?    // if this is true, ? then
      gradient_m*x + c_intercept  :                                       // zb[] equals this : else,
      x > length_until_runup - waterdepth/tan(theta1) ?                   // if this is true ? then
      tan(theta1) * (x-length_until_runup) :                              // zb[] equals this : else,
      -waterdepth;                                                        // zb[] equals this


    foreach_layer() {
      // Granular fluid
      if (x <= width/2) {
        // DEPTH
        h[] =  (height_cm)/(nl);
        // BUOYANCY
        T[] = rho_1;
      }

      // Granular fluid (sloped)
      else if (x < width) { 
        // DEPTH
        h[] = (-1*x*(height_cm/(width/2.)) + 2*height_cm)/nl;
        // BUOYANCY
        T[] = rho_1;
      }

      // Gap
      else if (x < width + slope_exposed_height/tan(theta0)) { 
        // DEPTH
        h[] = dry;  // "dry"
        // BUOYANCY
        T[] = rho_1;
      }

      // Water
      else {
        // DEPTH
        h[] = max(-zb[], dry)/nl;
        // BUOYANCY
        T[] = rho_0;
      }

      // Find max/min T for plotting.
      if(T[] > max_T) max_T = T[];
      if(T[] < min_T) min_T = T[];
    }
  }
}










/**
### Outputs */

event logfile (i += 10)
{
  static double rpe0 = 0., rpen = 0., tn = 0.;
  double rpe = RPE();
  double PE, KE;
  energy (&PE, &KE);
  if (i == 0) {
    rpe0 = rpe;
    fprintf (stderr, "t      dt      rpe-rpe0  rpe0    PE      KE        d_t(rpe)  mgp.i\n");
  }
  fprintf (stderr, "%-6g %-7g %-9.4g %-7.4g %-7.4g %-9.3g %-9.4g %-3d\n", t, dt,
	   rpe - rpe0, rpe0, PE, KE,
	   t > tn ? (rpe - rpen)/(t - tn)/L0 : 0.,
#if NH	   
	   mgp.i
#else
	   mgH.i
#endif
	   );
  rpen = rpe, tn = t;

  double max_runup_dist_tmp = 0.;
  foreach() {
    double H = 0.;
    foreach_layer() {
      H += h[];
    }
    if (H > runup_tolerance)
      if (x > max_runup_dist)
        max_runup_dist_tmp = x;
  }
  if (max_runup_dist_tmp > max_runup_dist)
    max_runup_dist = max_runup_dist_tmp;

}









/**
 * Output images including density, u.x...
*/

void setup (FILE * fp)
{
  fprintf (fp,
#if ISOPYCNAL
	   "set pm3d map corners2color c2\n"
#else
	   "set pm3d map\n"
#endif
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set xlabel 'x'\n"
	   "set ylabel 'depth'\n"
	   );
}






void plot_ux (FILE * fp)
{
  fprintf (fp,
	   "set cbrange [-0.1:5]\n" // u.x
	   "set title 't = %.2f'\n"
	   "sp '-' u 1:2:3\n", t);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}

void plot_density (FILE * fp)
{
  fprintf (fp,
	   "set cbrange [%.2f:%.2f]\n" // Density
	   "set title 't/T0 = %.2f'\n"
	   "sp '-' u 1:2:4\n", min_T, max_T, t);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}



event gnuplot (t += 0.1)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    setup (fp);
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1024,400\n"
	   "set xrange [0:12.5]\n"
	   "set yrange [%.3f:%.3f]\n"
	   "set size ratio -1\n"
	   "set ytics\n"
	   "set output 'gnuplot/closeup/plot-%06d.png'\n", -waterdepth-0.5, c_intercept+height_cm+0.5, i);
  plot_density (fp);
  
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1600,266\n"
	   "set xrange [0:%.3f]\n"
	   "set yrange [%.3f:%.3f]\n"
	   "set size ratio -1\n"
	   "set ytics\n"
	   "set output 'gnuplot/all/plot-%06d.png'\n", domain_length, -waterdepth-0.5, c_intercept+height_cm+0.5, i);
  plot_density (fp);

  fprintf (fp,
     "set term pngcairo font \",10\" size 1600,266\n"
     "set xrange [0:%.3f]\n"
     "set yrange [%.3f:%.3f]\n"
     "set size ratio -1\n"
     "set ytics\n"
     "set output 'gnuplot/all/ux/plot-%06d.png'\n", domain_length, -waterdepth-0.5, c_intercept+height_cm+0.5, i);
  plot_ux (fp);  
}






/**
 * Numerical gauges.
*/

Gauge gauges[] = {
  {"WG_10m",  10., 0.},
  {"WG_12m", 12., 0.},
  {"WG_14m", 14., 0.},
  {"WG_16m", 16., 0.},
  {"WG_18m", 18., 0.},
  {NULL}
};


event output (t += 0.1)
  output_gauges (gauges, {eta});



/**
 * At end, convert output images into animations.
*/

event end (t = 15)
{ 
  system ("for f in gnuplot/all/plot-*.png; do"
	  " convert $f ppm:- && rm -f $f; done | "
	  "ppm2mp4 movie_normal.mp4");
  system ("for f in gnuplot/all/ux/plot-*.png; do"
    " convert $f ppm:- && rm -f $f; done | "
    "ppm2mp4 movie_normal_ux.mp4");
  system ("for f in gnuplot/closeup/plot-*.png; do"
	  " convert $f ppm:- && rm -f $f; done | "
	  "ppm2mp4 movie_closeup.mp4");
  fprintf (stderr, "\n\nDone\n");
  fprintf (stderr, "\nMax Runup Distance:  %.4f\n", max_runup_dist);
}

