/**
# Capillary wave oscillation

Here, we test qualitatively the ability of the multilayer solver to
capture Laplace pressure and a capillary wave.  A periodic wave profil
is initialized. Surface is supposed to oscillate.

![Surface profile oscillation](oscillation/movie.mp4)

## Include */

#include "grid/multigrid1D.h"
#include "crobert/2_Implicit/hydro-tension.h"
//#include "layered/implicit.h"
#include "layered/nh.h"
#include "layered/remap.h"

/**
## Geometry and resolution*/
#define h_ 0.1
#define LEVEL 6
#define layers 1
#define T_END (10*2.*pi/sqrt(tanh(h_)))
const double h_0 = h_;
const double delta = 0.1;
const double k = 1.;

/**
## Main function
The test is done without any viscosity or gravity. */
int main()
{
  size (2.*pi);
  periodic (right);
  N = 1 << LEVEL;
  nl = layers;
  G = 0. [1,-2];
  NITERMAX = 1000;
  TOLERANCE = 1e-8 [*];
  linearised = true;
  system ("rm -f plot*.png");
  for (CFL_H = 10; CFL_H <= 200; CFL_H *= 2) {
    fprintf (stderr, "# CFL_H = %g\n", CFL_H);
    run();
    fprintf (stderr, "\n\n");
  }
}

/**
## Initialisation 
A wave relaxes slowly.*/
event init (i = 0)
{
  foreach() {
    double H = h_0*(1 + delta*cos(x*k));
    foreach_layer()
      h[] = H/nl;
  }
  boundary({h});
} 

/**
## Outputs*/
event outputs (i++)
{
  fprintf (stderr, "%g %g %g %g\n", t, statsf(eta).max,
	   h_0 + h_0 * delta * fabs(cos(1[0,-1]*t*sqrt(tanh(h_0*k)))),
	   dt/sqrt(cube(L0/N)/(4.*pi)));
}

/**
# Movie
To see the movement of the wave*/
event movie (i++)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio 0.3\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][%g:%g]'-' u 1:(-1):2 w filledcu lc 3 t ''",
	   i, t, X0, X0 + L0, (1. - delta)*h_0, (1. + delta)*h_0);
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

event moviemaker (t = T_END)
{
  system ("rm -f movie.mp4 && "
	  "ffmpeg -r 25 -f image2 -i plot%04d.png "
	  "-vcodec libx264 -vf format=yuv420p -movflags +faststart "
	  "movie.mp4 2> /dev/null && "
	  "rm -f plot*.png");
}

/**
Note that the theoretical solution is only valid in the linear limit,
which is not verified for the relative large amplitude we chose.

~~~gnuplot Height evolution
set xlabel 'time t'
set ylabel 'height h' 
set xrange [175:]
set key below
plot 'log' u 1:2 index 'CFL_H = 10' w lp t 'CFL_H = 10',	 \
     'log' u 1:2 index 'CFL_H = 20' w lp t 'CFL_H = 20',	 \
     'log' u 1:2 index 'CFL_H = 40' w lp t 'CFL_H = 40',	 \
     'log' u 1:2 index 'CFL_H = 80' w lp t 'CFL_H = 80',	 \
     'log' u 1:2 index 'CFL_H = 160' w lp t 'CFL_H = 160',	 \
     'log' u 1:3 index 'CFL_H = 10' w l t 'theory'
~~~
*/
