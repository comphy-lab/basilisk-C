/**
# Solitory wave surging over an emerged reef crest from Roeber and Cheung (2012) */

# include "grid/multigrid1D.h"
# include "green-naghdi.h"

/**
The domain is 102.4 metres long, slipt into 0.1m cells.
I tried a number of breaking values to turn off dispersion and found 0.4 worked best
for this case. gravity is set to 9.81 */

int main()
{
  L0 = 102.4;
  X0 = -18.7;
  G = 9.81;
  N = 1 << 10;
  breaking = 0.4;
  run();
}

/**
The initial conditions are for the Green-Naghdi soliton (soliton.c) in
a water depth of 2.5 metres and a relative soliton amplitude of 0.30. */

double h0 = 2.5, a = 0.3;

double sech2 (double x) {
  double a = 2./(exp(x) + exp(-x));
  return a*a;
}

double soliton (double x, double t)
{
  double c = sqrt(G*(1. + a)*h0), psi = x - c*t;
  double k = sqrt(3.*a*h0)/(2.*h0*sqrt(h0*(1. + a)));
  return a*h0*sech2 (k*psi);
}

event init (i = 0)
{
  double c = sqrt(G*(1. + a)*h0);
  foreach() {
    double eta = soliton (x + 5., t);
    
/** The bathymetry as outlined in Roeber and Cheung 2012. */
    x -= 0.;
    zb[] = (x < 25.9 ? -2.5 :
	    x < 56.65 ? min(-4.65 + (x *0.0833), 0.065) :
	    x < 57.95 ? 0.065 :
	    x < 60.9 ? max(4.87 + (x) *-0.083, -0.136) :
	    -0.136);
    h[] = max (0., eta-zb[]);
    u.x[] = c*eta/(h0 + eta);

  }
}

/**
Add friction and stop simulation at 50s. */

event timeseries (i++; t <= 50) {
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-2*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
  }
}

/**
gnuplot is used to to visualise the wave profile as the simulation
runs and to output a snapshot every 5 seconds $t=40$. */

void plot_profile (double t, FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "set xlabel 'x (m)'\n"
  	   "set ylabel 'y (m)'\n"
	   "p [0:83.7][-2.5:1]'-' u 1:3:2 w filledcu lc 3 t '',"
	   " '' u 1:(-2.5):3 t '' w filledcu lc 0\n", t);
  foreach()
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n\n");
  fflush (fp);
}
/** Un-comment below to visualise simulating during processing 

event gnuplot (t += 0.05) {
  static FILE * fp = popen ("gnuplot", "w");
  plot_profile (t, fp);
  fprintf (stderr, "%g %g\n", t, interpolate (eta, 17.3, 0.));
} */

/** To make a movie change to t += 0.05, then use: ffmpeg -f image2 -r 30 -i reefcrest/'snapshot-1%05d.png' ReefCrest.mp4

![Snapshot of plunge. The reef is white.](reefcrest/snapshot-100230.png)
*/

event gnuplot (t = 11.5) {
  FILE * fp = popen ("gnuplot", "w");
  fprintf (fp,
           "set term pngcairo enhanced size 640,200 font \",8\"\n"
           "set output 'snapshot-%.0f.png'\n", (t*20)+100000);
  plot_profile (t, fp);
  pclose (fp);
}

/** Un-comment to output profile data for pressure, surface and bathymetry every second

event profile (t += 1.00) {
  char name[20]; sprintf (name, "p-%g", t);
  FILE * fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g %g %g\n", x, h[], eta[], zb[]);
  fclose (fp);
}
*/

/**
Output profile data was used to compare basilisk results with results for Roeber and Cheung (2012)

![Basilisk Green-Naghdi solution (black line with blue fill) compared with the basilisk Saint Venant solution (grey line) and Roeber and Cheung (2012) solution (red line).](/ReefCrestOut.jpg)

Reference

Roeber, V. and Cheung, K.F., 2012. Boussinesq-type model for energetic breaking waves in fringing reef environments. Coastal Engineering, 70: 1-20. */

