/**
## MWE for surface relaxation due to an homogeneous surface tension
A periodic wave profil is initialized. Surface tension is supposed to flatten the surface.

![Surface profil](relaxation/movie.mp4)
*/

double initial_mass = 0.;
#include "grid/multigrid1D.h"
#include "../hydro_c.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "../solute.h"
#include "../surface.h"

/**
## Geometry and resolution*/
#define L 10.
#define h_ 0.5
double ak = 0.35;
#define LEVEL 7
#define layers 10
#define T_END 20.
#define DELTA_T (T_END/100.)
#define MY_TOLERANCE 1e-11

/**
## Constants*/
#define k_  (2.*pi)

/**
## Physical parameters*/
#define g_ 0.
#define Re 1.

/**
We set the value of the contact angle of the liquid on the walls. */
#define liquid_theta (90.*pi/180.)

/**
## Main function
There is no diffusion, no adsorption, no viscosity, no gravity*/
int main()
{
  L0 = L;
  origin (0.);
  periodic(right);  
  /**theta_c = liquid_theta;*/
  N = 1 << LEVEL;
  TOLERANCE = MY_TOLERANCE;
  nl = layers;
  G = 0;
  nu = 1/Re;
  run();
}

/**
## Boundary conditions

We set the boundary conditions on the height functions thanks to contact.h. */
/**H[left] = contact_angle (liquid_theta, L0/N);
H[right] = contact_angle (liquid_theta, L0/N);
*/
/**
## Initialisation 
A wave is convected from left to right at a constant velocity. Bulk and surface concentrations are inhomogeneous*/
#include "test/stokes.h"

event init (i = 0)
{
  scalar c, h;
  for(c, h in cl, hl)
    foreach(){
      double H = h_ + 4*wave(x/2.5, 0);
      c[] = 0;
      h[] = H/nl;
      gam[] = 1.;
    }
  heights(hl, H);
}

/**
## Outputs*/
event outputs(t += DELTA_T; t = 0.; t <= T_END) {
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