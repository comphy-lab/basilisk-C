/**
# Meniscus
In this example we examine the ability of the explicit hydro-capillary solver to recover equilibrium interface shapes. We initialize the interface shape with a circular arc and check that it remains so through time.
*/

#include "grid/multigrid1D.h"
#include "./hydrocapillary_explicit.h"
#include "layered/remap.h"

/**
## Wall boundary conditions

At the wall we impose a contact angle and also wish to make sure that the Laplace pressure flux is zero. This can be achieved thanks to the second row of ghost cells (which can be accessed with [Antoon's trick](/sandbox/Antoonvh/two_layers.c)).
*/
#define BGHOSTS 2

@define ghost_layer (neighbor.i < GHOSTS ? (GHOSTS - neighbor.i) : neighbor.i - (1 << level) - GHOSTS + 1) // gives 1 for the first layer, 2 for the second, deepest, layer
@define compass (neighbor.i < GHOSTS ? 1 : -1) // a compass to distinguish between left from right
@define slopebnd ((val(_s,0,0,0)-val(_s,compass,0,0))/Delta)
//@define slopebnd(th) (1./tan(th))
@define tgtdiff(th) (2.*cos(th)-slopebnd/sqrt(1+sq(slopebnd)))
//@define wall(th) (val(_s,0,0,0) + Delta / tan(th) + (ghost_layer - 1)*Delta*tgtdiff(th)/clamp(sqrt(1-sq(tgtdiff(th))),1,1))
  @define wall(th) (val(_s,0,0,0) + Delta / tan(th) + (ghost_layer - 1)*Delta*tgtdiff(th))
//  @define wall(th) (val(_s,0,0,0) + Delta / tan(th) + (ghost_layer - 1)*Delta*tgtdiff(th)/sin(th))
//  @define wall(th) (val(_s,0,0,0) + (ghost_layer)*Delta*slopebnd)
//  @define wall(th) (val(_s,0,0,0) + Delta / tan(th) + (ghost_layer - 1)*Delta / tan(th) )

/**
## Geometry
*/

#define L 1.
#define Tend 1.
#define DELTA_T (Tend/400.)
#define LEVEL 7
#define layers 4
#define contact_angle (60.*pi/180.)

/**
We define the initial circular shape of the interface with
*/

#define radius (L/(2.*cos(contact_angle)))
#define zc (radius*sin(contact_angle))
#define meniscus_shape(x) (1. + zc-pow(sq(radius)-sq(x-L/2.),0.5))
//#define meniscus_shape(x) (1. + 0.2*cos(2*pi*x))

/**
## Main function
We run this test without gravity nor viscosity.
*/

int main() {
  CFL_H = 0.5;
  L0 = L;
  origin (0.);
  G = 0.;
  nl = layers;
  gradient = NULL; // this allows for "steep" slopes
  N = 1 << LEVEL;
  run();
}

/**
## Initialization
We make use of the previously defined `wall` boundary condition.
*/

event init (i = 0) {
  eta[left]  = neumann (1./tan(contact_angle));
  eta[right] = neumann (1./tan(contact_angle));
  curvature[left] = neumann(0);
  curvature[right] = neumann(0);  
  foreach() {
    double H = meniscus_shape(x);
    foreach_layer() {
      h[] = H/nl;
      sigma[] = 1.;
    }
  }
  boundary ({eta, sigma});
}

/**
## Movie
*/

event outputs(t = 0.; t += DELTA_T; t <= Tend) {
  static FILE * fpm = fopen ("height.txt", "w");
  fprintf (fpm, "%g %g %g\n", t, statsf(eta).max, 1.);
}

event logfile (i++; t <= Tend)
{
}

#if 1
event movie (t = 0.; t += DELTA_T; t <= Tend)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio -1\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][0:2]'-' u 1:(-1):2 w filledcu lc 3 t ''",
	   i/2, t, X0, X0 + L0);
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
#endif
