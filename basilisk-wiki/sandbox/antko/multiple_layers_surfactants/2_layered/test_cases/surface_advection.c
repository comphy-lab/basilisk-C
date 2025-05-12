/**
## MWE for surface convection
A wave profil in a domain with periodic conditions is convected with a fix velocity from left to right.
The convection of surfacic concentration is monitored.
*/

double initial_mass = 0.;
#include "grid/multigrid1D.h"
#include "../hydro_c.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "../surface.h"
#include "../solute.h"

/**
## Geometry and resolution*/
#define L 1.
#define h_ 0.5
double ak = 0.35;
#define LEVEL 7
#define layers 10
#define T_END 15.
#define DELTA_T (T_END/1000.)
#define MY_TOLERANCE 1e-11

/**
## Constants*/
#define k_  (2.*pi)

/**
## Physical parameters*/
#define g_ 0.
#define Re 1.
#define solute0 1.
#define surface0 1.

/**
## Main function
There is no diffusion, no adsorption, no viscosity, no gravity*/
int main()
{
  L0 = L;
  origin (0.);
  periodic (right);
  N = 1 << LEVEL;
  TOLERANCE = MY_TOLERANCE;
  nl = layers;
  G = g_;
  run();
}


/**
## Initialisation 
A wave is convected from left to right at a constant velocity. Bulk and surface concentrations are inhomogeneous*/
#include "test/stokes.h"

event init (i = 0)
{
  scalar h, c;
  vector u;
  for (u, h,c in ul, hl,cl) {
    foreach() {
      double H = h_ + wave(x, 0);
      h[] = H/nl;
  	  c[] = solute0 * (1 + wave(x, 0));
  	  u.x[] = 1.;
    }
    boundary({c});
  }
  interface_area(hl, area);
  foreach()
    surface[] = surface0*area[] * (1 + wave(x, 0));
}

/**
## Outputs*/
event outputs(t += DELTA_T; t = 0.; t <= T_END) {
  static FILE * fpSurf = fopen("dataSurf.txt", "w");
  static FILE * fp = fopen("surf_mass.txt", "w");
  
  double total_solute = 0 ;
  double total_surface = 0 ;
  double min_solute = nodata;
  double max_solute = 0;
  double min_surface = nodata;
  double max_surface = 0;
  
  scalar h, c;
  foreach (){
    for (h,c in hl,cl) {
      total_solute += h[]*c[]*Delta;
      min_solute = min(min_solute, h[]*c[]);
      max_solute = max(max_solute, h[]*c[]);
    }
    total_surface += surface[]*Delta;
    min_surface = min(min_solute, surface[]);
    max_surface = max(max_solute, surface[]);    
    fprintf (fpSurf, "%g %g %g %g %g \n", t, x, surface[], area[], surface[]/area[]);
  }
  initial_mass = (initial_mass == 0) ? total_solute + total_surface : initial_mass;
  fprintf (fp, "%g %g %g %g %g %g %g %g\n", t, total_solute, total_surface, total_solute + total_surface - initial_mass, min_solute, max_solute, min_surface, max_surface);
}

/**
# Movie
To see the movement of the wave
event movie (i += 3)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio -1\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][0:1]'-' u 1:(-1):2 w filledcu lc 3 t ''",
	   i/3, t/(k_/sqrt(k_)), X0, X0 + L0);
  fprintf (fp, "\n");
  foreach_leaf() {
    double H = 0.;
    for (scalar h in hl)
      H += h[];
    fprintf (fp, "%g %g %g", x, zb[] + H, zb[]);
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}

event moviemaker (t = end)
{
  system ("rm -f movie.mp4 && "
	  "ffmpeg -r 25 -f image2 -i plot%04d.png "
	  "-vcodec libx264 -vf format=yuv420p -movflags +faststart "
	  "movie.mp4 2> /dev/null && "
	  "rm -f plot*.png");
}
*/

/** ##Results
~~~gnuplot Difference between total mass and initial total mass
set terminal @PNG enhanced size 320,320 font ",8"
set output 'sum.png'
set xlabel "t"
set key left
plot \
  './surf_mass.txt' u 1:4 t 'sum' w l
~~~
Total mass is well conserved, the error is the machine error.

~~~gnuplot Minimum surface concentration.
set output 'min.png'
set xlabel "t"
set key left
plot \
  './surf_mass.txt' u 1:5 t 'sum' w l
~~~

~~~gnuplot Maximum surface concentration.
set output 'max.png'
set xlabel "t"
set key left
plot \
  './surf_mass.txt' u 1:6 t 'sum' w l
~~~
There is some unwanted diffusion for minimum and maximum surface. This error reduces when Delta increases. A better surface area approximation (not linear) could help to avoid this diffusion.
*/
