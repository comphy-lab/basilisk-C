/**
# Capillary wave oscillation
**Fork from crobert sandbox**
Here, we test qualitatively the ability of the multilayer solver to capture Laplce pressure and a capillary wave.
A periodic wave profil is initialized. Surface is supposed to oscillate.

![Surface profil oscillation](oscillation/movie.mp4)
*/
 
/**
## Include*/
#include "grid/multigrid1D.h"
#include "./hydrocapillary.h"
#include "layered/remap.h"
#include "./hydrocapillary_implicit.h"
//#include "../nh.h"

/**
## Geometry and resolution*/
#define h_ 0.1
#define ah h_*0.01
#define LEVEL 6
#define layers 1
#define T_END (10*2.*pi/sqrt(tanh(h_)))
#define DELTA_T (T_END/400.)

/**
## Main function
The test is done without any viscosity or gravity. */
int main()
{
  CFL_H = 5.;
  size (2.*pi);
  periodic (right);
  N = 1 << LEVEL;
  nl = layers;
  G = 0.;
  TOLERANCE = 1e-8;
  run();
}

/**
## Initialisation 
A wave relaxe slowly.*/
event init (i = 0)
{
  foreach() {
    double H = h_ + ah*cos(x);
    foreach_layer() {
      h[] = H/nl;
      sigma[] = 1.;
    }
  }
  boundary({h, sigma});
} 

/**
## Outputs*/
event outputs(t = 0.; t += DELTA_T; t <= T_END) {
  static FILE * fpm = fopen ("height.txt", "w");
  fprintf (fpm, "%g %g %g\n", t, statsf(eta).max, h_ + ah * fabs(cos(t*sqrt(tanh(h_)))));
}


/**# Movie
To see the movement of the wave*/
event movie (t = 0.; t += DELTA_T; t <= T_END)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio 0.3\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][0.099:0.101]'-' u 1:(-1):2 w filledcu lc 3 t ''",
	   (int)(t/DELTA_T), t, X0, X0 + L0);
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
     'height.txt' u 1:3 w l t 'model' lt 2 lw 2
~~~
*/
