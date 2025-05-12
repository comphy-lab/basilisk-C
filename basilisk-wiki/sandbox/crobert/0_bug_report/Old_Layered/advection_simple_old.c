/**
## MWE for surface convection
A wave profil in a domain with periodic conditions is convected with a constant velocity from left to right. The convection of surfacic concentration and the height are monitored. The maximum values should remain constant.

![Animation of the free-surface](advection_simple_old/movie.mp4)
*/

#include "grid/multigrid1D.h"
#include "./hydro_c.h"
#include "./nh_c.h"
#include "./remap.h"

/** 
## Geometry and resolution*/
#define L 8.
#define h_ 0.5
double ak = 0.35;
#define LEVEL 8  
#define layers 10
#define T_END 10.
#define DELTA_T (T_END/100.)
#define MY_TOLERANCE 1e-11
 
/**
## Constants*/
#define k_  (2.*pi)

/**
## Physical parameters*/
#define g_   0.


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
  
  /**
  The default minmod slope limiter cuts off wave crests. This can be fixed
  by turning off limiting. */ 
  
  gradient = NULL;
  run();
}

/**
## Initialisation 
A wave is convected from left to right at a constant velocity. Bulk and surface concentrations are inhomogeneous*/
#include "test/stokes.h"

event init (i = 0)
{
  scalar h;
  vector u;
  for (u, h in ul, hl) {
    foreach() {
      double H = h_ + ak/5.*sin(2*pi*x);
        h[] = H/nl;
        u.x[] = 1.;
    }
  }
}

/**
## Outputs*/
event outputs(t += DELTA_T; t = 0.; t <= T_END) {
}

event outputsTime(i++) {
  fprintf(stdout, "i = %i ; t = %g\n", i, t);
}

/**
# Movie
To see the movement of the wave */
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