/**
# Kapitza waves
 
![Surface Profil](Kapitza/movie.mp4)
*/

/**
## Include and parameters
This test uses the multilayer solver and the surface tension additive. */
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"

#define beta (pi/180 * 6.4)

#define h_ (1./cos(beta))
#define Q_ 1.

#define Re 19.33
#define Fr2 (Re*sin(beta)/3.)
#define We 5.43
#define F 0

#define L 500
#define T_END 500.
#define DELTA_T (T_END/500.)

/**
## Main function
The boundary conditions are set to periodic. The test is done with zero-gravity. */
int main()
{
  L0 = L;
  N = 1024;
  nl = 20;
  
  G = 1./Fr2;
  nu = 1./Re;
  
  run();
}

/**
## Initialisation */

/**
The inflow is a parabolic profile with a total flow rate Q. The function below computes the height zc of the middle of the layer and returns the corresponding velocity. For the multilayer solver, we need some trickery to define the inflow velocity profile as a function of zz. We need to sum the layer thicknesses to get the total depth H and the height zc of the current layer (of index point.l).*/
double Qparab (Point point, double Q)
{
  double H = 0.;
  double zc = 0.;
  for (int l = 0; l < nl; l++) {
    H += h[0,0,l];
    if (l < point.l)
      zc += h[0,0,l];
  }
  zc += (h[]/2.);
  return 3./2.*Q*(2.*(zc/H) - sq(zc/H));
}

double A(double t)
{
return (1 + 0.05*sin(2*pi*F*t));
}

event init (i = 0)
{
  foreach() {
    zb[] = (L - x) * sin(beta);
    sigma[] = We;
    double H = h_;
    foreach_layer()
      h[] = H/nl;
    foreach_layer()
      u.x[] = Qparab (point, Q_)/h_;
  }
  
  u.n[left] = Qparab (point, Q_)/h_;
  hu.n[left] = Qparab (point, Q_)/nl;
  h[left] = h_/nl;
  eta[left] = ((L + Delta/2)*sin(beta) + h_);

  u.n[right] = neumann(0.);
  eta[right] = (- Delta/2 *sin(beta) + h_*0.9);
  phi[right] = dirichlet(0.);
  boundary(all);
}

/**
## Outputs*/
event outputs(t += DELTA_T; t = 0.; t <= T_END) {
}

/**
# Movie 
This is how we export the headline movie. */
event movie (i +=  100)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio 0.2\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][0.5:1.8]'-' u 1:3:2 w filledcu lc 3 t ''",
	   i/100, t, X0, X0 + L0);
  fprintf (fp, "\n");
  foreach_leaf() {
    double H = 0.;
    foreach_layer()
      H += h[];
    fprintf (fp, "%g %g 0", x, H);
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
     

  static FILE * fpp = fopen ("note", "w");
  foreach() {
    double H = 0;
    double Q = 0;
    foreach_layer() {
      H += h[];
      Q += h[]*u.x[];
    }
    fprintf (fpp, "%g %d %g %g\n", x, point.l, H, Q);  
  }
}
