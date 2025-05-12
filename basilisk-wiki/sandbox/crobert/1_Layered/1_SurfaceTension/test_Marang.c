/**
# Marangoni flow vs Laplace Pressure

We test here the ability of the multilayer solver to reproduce correctly the Marangoni flow given a surface tension gradient. 
This is a qualitative test controlled by surface tension and under zero-gravity.

A periodic tank of finite dimensions containing a liquid film of mean height $h_0$ is considered. 
A periodic surface tension is imposed on the surface of the liquid (phisically, such a gradient could be obtained for a thin film with a temperature gradient).
Marangoni flows tend to accumulate liquid on the low surface tension areas whereas Laplace pressure tends to flatten the surface.
 
![Surface Profil](test_Marang/movie.mp4)
*/

/**
## Include and parameters
This test uses the multilayer solver and the surface tension additive. */
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"

#define h_ 0.5
#define T_END 5.
#define DELTA_T (T_END/100.)

/**
## Main function
The boundary conditions are set to periodic. The test is done with zero-gravity. */
int main()
{
  periodic(right);
  N = 32;
  nl = 2;
  G = 0;
  nu = 10;
  run();
}

/**
## Initialisation 
A wave relaxe slowly.*/
event init (i = 0)
{
  foreach() {
    double H = h_;
    foreach_layer() {
      h[] = H/nl;
      sigma[] = 1. + 0.5*cos(2*pi*x);
    }
  }
}

/**
## Outputs*/
event outputs(t += DELTA_T; t = 0.; t <= T_END) {
}

/**
# Movie 
This is how we export the headline movie. */
event movie (i += 10)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio .4\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
     "set xlabel 'x'\n"
     "set ylabel 'height'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][0.2:0.6]'-' u 1:(-1):2 w filledcu lc 3 t ''",
	   i/10, t, X0, X0 + L0);
  fprintf (fp, "\n");
  foreach_leaf() {
    double H = 0.;
    foreach_layer()
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