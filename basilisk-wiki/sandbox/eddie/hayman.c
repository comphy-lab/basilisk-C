/**
# Wave transformation across Hayman Island reef, from Gourlay (1994) */

# include "grid/multigrid1D.h"
# include "green-naghdi.h"
/** change DEPTH (inner reef flat), HT (wave height) TS (wave period) for each scenario */
#define DEPTH 1.4
#define HS 3.03
#define TS 5.8

int main()
{
  L0 = 1024;
  G = 9.81;
  N = 1 << 11;
  breaking = 1;
/** to minamise dissipation turn off limiting
  gradient = NULL; */
  run();
}

scalar hmax[];
scalar SH[];
scalar SHc[];
scalar AH[];
scalar Vel[];
scalar SVel[];
scalar AVel[];
scalar vmax[];
scalar vmin[];

u.n[left]  = - radiation ((HS)*sin(2.*pi*t/TS));
u.n[right] = + radiation (0);

event init (i = 0)
{
/** define Hayman reef bathymetry as described in Gourlay (1994)*/
foreach() {
    x += 0.;
    zb[] = (
	    x < 253 ? max (min (-20 + ((x-180) *(1/4.5)), -3.8),-20) :
	    x < 276 ? min (-3.8 + ((x-256) *(0.03)), -3.2) :
	    x < 435 ? min (-3.2 + ((x-275) *(0.007)), -2.1) :
	    -2.1) + 2.1 - DEPTH;

    h[] = max (0., -zb[]);
    eta[] = h[] > dry ? h[] + zb[] : 0;
    boundary ({eta});
  }
}
/** Add implicit quadratic bottom firction, with extra friction before right boundary to eliminate reflection*/
event friction (i++) {
  foreach() 
       {
      double a = x < 512 ? h[] < dry ? HUGE : 1. + 0.01*dt*norm(u)/h[] :
		           h[] < dry ? HUGE : 1. + 0.08*(x - 512.)*dt*norm(u)/h[];
      foreach_dimension()
	u.x[] /= a;
	Vel[] = norm(u);
 }
}

/** calculate time averaged water level and velocity stats after wave field saturates reef*/
event timeseries (t = 200; i++) {
  foreach() {
 /** Set hmax as highest wave amplitude */
    if (h[] > dry && h[] + zb[] > hmax[])
      hmax[] = h[] + zb[];
 /** Set SH as the sum water level */
    if (h[] > dry)
      SH[] = SH[] + (h[] + zb[]);
 /** Set SHc as the number of sum water level */
    if (h[] > dry)
      SHc[] = SHc[] + 1;
 /** Set AH as the average water level */
    if (h[] > dry)
      AH[] = SH[] / SHc[];
 /** Set SVel as the sum velocity */
    if (h[] > dry)
      SVel[] = SVel[] + (Vel[]);
 /** Set AVel as the average velocity */
    if (h[] > dry)
      AVel[] = SVel[] / SHc[];
/** Set vmax as highest velocity */
    if (h[] > dry && Vel[] > vmax[])
      vmax[] = Vel[];
/** Set vmax as highest velocity */
    if (h[] > dry && Vel[] < vmin[])
      vmin[] = Vel[];
  }
}

void plot_profile (double t, FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "set xlabel 'x (m)'\n"
  	   "set ylabel 'y (m)'\n"
	   "p [200:512][-5:3]'-' u 1:3:2 w filledcu lc 3 t '',"
	   " '' u 1:(-5):3 t '' w filledcu lc -1\n", t);
  foreach()
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n\n");
  fflush (fp);
}
/** Visualise outputs during simulation 

event gnuplot (t += 0.1) {
  static FILE * fp = popen ("gnuplot", "w");
  plot_profile (t, fp);
  fprintf (stderr, "%g %g\n", t, interpolate (eta, 17.3, 0.));
} */

/** Save snapshot at end. To make a movie use t += 0.1, then: ffmpeg -f image2 -r 30 -i hayman/'snapshot-1%04d.png' hayman.mp4 
 
![Snapshot of waves. The reef is white.](hayman/snapshot-11500.png)
*/

event gnuplot (t = 300) {
  FILE * fp = popen ("gnuplot", "w");
  fprintf (fp,
           "set term pngcairo enhanced size 640,200 font \",8\"\n"
           "set output 'snapshot-%.0f.png'\n", (t*5)+10000);
  plot_profile (t, fp);
  pclose (fp);
} 

/** extract time series water level data every 20 meters across reef */
Gauge gauges[] = {
 {"g20", 20},
 {"g40", 40},
 {"g60", 60},
 {"g80", 80},
 {"g100", 100},
 {"g120", 120},
 {"g140", 140},
 {"g160", 160},
 {"g180", 180},
 {"g200", 200},
 {"g220", 220},
 {"g240", 240},
 {"g260", 260},
 {"g280", 280},
 {"g300", 300},
 {"g320", 320},
 {"g340", 340},
 {"g360", 360},
 {"g380", 380},
 {"g400", 400},
 {"g420", 420},
 {"g440", 440},
 {"g460", 460},
 {"g480", 480},
 {"g500", 500},
 {"gourlay250", 250},
 {"gourlay260", 260},
 {"gourlay270", 270},
 {"gourlay290", 290},
 {"gourlay370", 370},
 {"gourlay450", 450},
 {NULL}
};

event output (t += 0.1; t <= 300)
  output_gauges (gauges, {eta});

/** Output profile data for water level, setup, bathy and velocity */
event profile (t += 100) {
  char name[20]; sprintf (name, "p-%g", t);
  FILE * fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g %g %g %g %g %g %g %g %g\n", x, h[], u.x[], zb[], hmax[], SH[], AH[], eta[], Vel[], AVel[]);
  fclose (fp);
}
/** Reference
Gourlay MR (1994) Wave Transformation on a Coral-Reef, Coastal Engineering. 23:17-42.
Link to animation
https://vimeo.com/122389966
*/