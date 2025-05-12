/**
# Relaxation due to Laplace pressure

Here, we test qualitatively the ability of the multilayer solver to
capture the Laplace pressure. A periodic wave profil is
initialized. Surface tension is supposed to flatten the surface. A
very large timestep is used to check the stability of the implicit
scheme.

![Surface profil relaxation](relaxation/movie.mp4)
*/
 
/**
## Include*/
#include "grid/multigrid1D.h"
#include "crobert/2_Implicit/hydro-tension.h"
#include "layered/implicit.h"
#include "layered/remap.h"
#include "layered/perfs.h"

/**
## Geometry and resolution*/
#define L 1.
#define h_ 0.5
#define ah (h_/2.)
#define LEVEL 8
#define layers 3
#define T_END 3.
const double h_0 = h_ * 1.[1];
const double delta = 0.5;

/**
## Constants*/
const double k = 2.*pi;

/**
## Physical parameters*/
const double g = 0.;
const double Re =  0.01;

/**
## Main function
The test is done without any viscosity or gravity. */
int main()
{
  L0 = L * 1.[1];
  periodic (right);
  TOLERANCE = 1e-6;
  CFL_H = 100.;
  N = 1 << LEVEL;
  nl = layers;
  G = g;
  nu = 1/Re;
  system ("rm -f plot*.png");
  run();
}

/**
## Initialisation 
A wave relaxes slowly. */

event init (i = 0)
{
  foreach() {
    double H = h_0 * (1 + delta*cos(k*x));
    foreach_layer()
      h[] = H/nl;
  }
}

/**
## Outputs*/
event outputs (t <= T_END; i++) { 
  static FILE * fpm = fopen ("height.txt", "w");
  fprintf (fpm, "%g %g %g\n", t, statsf(eta).max, 0.5 + 0.1 * exp(- t*0.6 [0,-1]));
}

/**# Movie
To see the movement of the wave*/
event movie (i += 10)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio -1\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][0:1]'-' u 1:(-1):2 w filledcu lc 3 t ''",
	   i/10, t, X0, X0 + L0);
  fprintf (fp, "\n");
  foreach (serial) {
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

/**
~~~gnuplot Height evolution
set terminal svg enhanced size 640,640 font ",8"
set xlabel 'time t'
set ylabel 'height h' 
plot 'height.txt' u 1:2 w p t 'maximum height' pt 3 ps 1.0, \
     0.5 + 0.1*exp(-0.66*x) t 'model' lt 2 lw 4
~~~
*/
