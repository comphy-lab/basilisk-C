/**
# Tsunami generation by a "pyroclastic" density current

This is a simulation of a pyroclastic density current entering water
and generating tsunami waves. There is also an option for a dense
layer at base of flow, overlayed by a lighter layer (not shown in this
version). The setup is dimensionless and based on the experimental
geometry of [Bougouin et al., 2020](#bougouin2020), where a laboratory
fluidised granular flow runs down a ramp and into water. This is the
multilayer version of the Navier-Stokes case ([Battershill et al.,
2021](#battershill2021)).

The setup is adapted from [overflow.c](/src/test/overflow.c) to
simulate a flow of variable density layers into water, using the
[Boussinesq buoyancy](/src/layered/dr.h) module of the
[non-hydrostatic](/src/layered/nh.h) [layered
scheme](/src/layered/hydro.h).

![Closeup view: the colorscale is the relative density
 variation.](slide/movie_closeup.mp4)(autoplay loop width="100%")

![Full tank view: the colorscale is the relative density
 variation.](slide/movie_normal.mp4)(autoplay loop width="100%")

*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"

/**
Setting $\Delta \rho$ directly.
*/

#define drho(T) T

/**
Including the relevant modules and intialising the simulation
*/

#include "layered/dr.h"

#include "layered/remap.h"
#include "layered/rpe.h"
#include "layered/perfs.h"

#define nlayers 24

double rho_0 = 0.01;
double rho_2 = -0.01;
double rho_3 = 0.;

//double nu_H = 0.0001; // 1e-2;

int main()
{
  size (32);
  //breaking = 0.7;
  N = 512; //2048
  nl = nlayers;
  nu = 0.0001;
  DT = 0.0025;
  G = 1.;
  cell_lim = mono_limit;
  system ("rm -r gnuplot");
  system ("mkdir gnuplot");
  system ("mkdir gnuplot/all");
  system ("mkdir gnuplot/closeup");
  system ("mkdir gnuplot/closeup/density");
  system ("mkdir gnuplot/closeup/ux");
  const scalar slip[] = HUGE;
  lambda_b = slip;
  run();
}

/**
Define the geometry. 
*/

#define theta0 0.268
#define theta1 0.55
#define width 0.338/0.265
#define gradient_m -1.*tan(theta0)
#define Hi 1.
#define slope_exposed 1.
#define Hel (slope_exposed/0.265)*sin(theta0)
#define c_intercept tan(theta0)*(width+ (Hi+Hel)/tan(theta0)) - Hi
#define c_intercept2 Hel + (1./tan(theta1))*width 
#define height_cm 0.185/0.265

/**
Initialise domain.
*/

event init (i = 0)
{
  foreach() {
    zb[] =  x < width + (Hi+Hel)/tan(theta0) ?
      gradient_m*x + c_intercept  : -Hi; 
    foreach_layer() {
      if (x <= width/2) {// Granular fluid
	// DEPTH
	h[] =  (height_cm)/(nl);
	// BUOYANCY
	T[] = rho_0;
      }
      else if (x < width) { //Granular fluid
	h[] = (-1*x*(height_cm/(width/2.)) + 2*height_cm)/nl;
	T[] = rho_0;
      }
      
      else if (x < width + Hel/tan(theta0)) { // Gap
	// DEPTH
	h[] = 0.00002;
	// BUOYANCY
	T[] = rho_0;
      }
      else { //Water
	// DEPTH
	h[] = - zb[]/nl;
	// BUOYANCY
	T[] = rho_3;
      }
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
    fprintf (stderr, "t  dt   rpe-rpe0   rpe0   PE   KE   d_t(rpe)   mgp.i\n");
  }
  fprintf (stderr, "%g %g %.12g %.12g %.12g %.12g %.12g %d\n", t, dt,
	   rpe - rpe0, rpe0, PE, KE,
	   t > tn ? (rpe - rpen)/(t - tn)/L0 : 0.,

#if NH	   
	   mgp.i
#else
	   mgH.i
#endif
	   );
  rpen = rpe, tn = t;
}

/**
 Animations including density, u.x...
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
	   "set cbrange [-0.1:30]\n" // u.x
	   "set title 't = %.2f'\n"
	   "sp '-' u 1:2:3\n",t);
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
	   "set cbrange [-0.01:0.01]\n" // Density
	   "set title 't/T0 = %.2f'\n"
	   "sp '-' u 1:2:4\n",t);
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
  // Uncomment for animation during run.
  /*
    if (getenv ("DISPLAY")) {
    fprintf (fp, "set term x11\n");
    plot_density (fp);
    }
  */
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1024,400\n"
	   "set xrange [0:12.5]\n"
	   "set yrange [-1.5:2.5]\n"
	   "set size ratio -1\n"
	   "set ytics\n"
	   "set output 'gnuplot/closeup/density/plot-%06d.png'\n", i);
  plot_density (fp);
  
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1600,266\n"
	   "set xrange [0:32]\n"
	   "set yrange [-1.5:2.5]\n"
	   "set size ratio -1\n"
	   "unset ytics\n"
	   "set output 'gnuplot/all/plot-%06d.png'\n", i);
  plot_density (fp);  
}

/**
Equally placed wave gauges.
*/

Gauge gauges[] = {
  {"WG_3m",  3./0.265, 0.},
  {"WG_4m", 4./0.265, 0.},
  {"WG_5m", 5./0.265, 0.},
  {"WG_6m", 6./0.265, 0.},
  {"WG_7m", 7./0.265, 0.},
  {NULL}
};

event output (t += 0.1)
  output_gauges (gauges, {eta});

/**
We plot the wave gauges:

~~~gnuplot Evolution of the free surface elevation at x = 
set xlabel 't/T0'
set ylabel 'y'
plot 'WG_3m' u 1:2 w l t 'x = 11.3','WG_4m' u 1:2 w l t 'x = 15.1', 'WG_5m' u 1:2 w l t 'x = 18.9'
~~~
*/

event end (t = 25.)
{ 
  system ("for f in gnuplot/all/plot-*.png; do"
	  " convert $f ppm:- && rm -f $f; done | "
	  "ppm2mp4 movie_normal.mp4");
  system ("for f in gnuplot/closeup/density/plot-*.png; do"
	  " convert $f ppm:- && rm -f $f; done | "
	  "ppm2mp4 movie_closeup.mp4");
  fprintf (stderr, "\n\nDone\n");
}

/**
## References

~~~bib
@article{bougouin2020,
  title={Impact of fluidized granular flows into water: implications for tsunamis generated by pyroclastic flows},
  author={Bougouin, Alexis and Paris, Raphael and Roche, Olivier},
  journal={Journal of Geophysical Research: Solid Earth},
  volume={125},
  number={5},
  pages={e2019JB018954},
  year={2020},
  publisher={Wiley Online Library},
  doi={10.1029/2019JB018954}
}

@article{battershill2021,
  title={Numerical simulations of a fluidized granular flow entry into water: insights into modeling tsunami generation by pyroclastic density currents},
  author={Battershill, Lily and Whittaker, Colin and Lane, Emily and Popinet, Stephane and White, James and Power, William and Nomikou, P},
  year={2021}
}
~~~
*/

