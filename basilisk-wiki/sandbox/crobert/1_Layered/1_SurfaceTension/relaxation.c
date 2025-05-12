/**
# Relaxation due to Laplace pressure
Here, we test qualitatively the ability of the multilayer solver to capture Laplce pressure.
A periodic wave profil is initialized. Surface tension is supposed to flatten the surface. 

![Surface profil relaxation](relaxation/movie.mp4)
*/
 
/**
## Include*/
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"

/**
## Geometry and resolution*/
#define L 1.
#define h_ 0.5
#define ah h_/5
#define LEVEL 8
#define layers 3
#define T_END 10.
#define DELTA_T (T_END/100.)

/**
## Constants*/
#define k_  (2.*pi)

/**
## Physical parameters*/
#define g_ 0.
#define Re 0.01

/**
## Main function
The test is done without any viscosity or gravity. */
int main()
{
  L0 = L;
  periodic(right);
  N = 1 << LEVEL;
  nl = layers;
  G = g_;
  nu = 1/Re;
  run();
}

/**
## Initialisation 
A wave relaxe slowly.*/
event init (i = 0)
{
  foreach() {
    double H = h_ + ah*cos(k_*x);
    foreach_layer() {
      h[] = H/nl;
      sigma[] = 1.;
    }
  }
}

/**
## Outputs*/
event outputs(t = 0.; t += DELTA_T; t <= T_END) {
  if (i < 10) {
    static FILE * fpp = fopen ("phi.txt", "w");
    foreach() {
      fprintf (fpp, "phi_th : %.4g", phi_s[]);
      foreach_layer ()
        fprintf (fpp, "      %.4g       ", phi[]);
      fprintf (fpp, "\n");
    }
    fprintf (fpp, "\n");
  }
  
  static FILE * fpm = fopen ("height.txt", "w");
  fprintf (fpm, "%g %g %g\n", t, statsf(eta).max, 0.5 + 0.1 * exp(- t*0.6));
}


/**# Movie
To see the movement of the wave*/
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
	   i/3, t, X0, X0 + L0);
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

/**
~~~gnuplot Height evolution

set terminal svg enhanced size 640,640 font ",8"
set xlabel 'time t'
set ylabel 'height h' 
plot 'height.txt' u 1:2 w p t 'maximum height' pt 3 ps 1.0, \
     'height.txt' u 1:3 w l t 'model' lt 2 lw 4
~~~
*/