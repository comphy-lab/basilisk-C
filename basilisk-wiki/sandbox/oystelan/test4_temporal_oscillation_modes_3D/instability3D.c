/**
# Water at rest over a relatively steep sloping bottom (on the diagonal) (3D)

This testcase triggers temporal oscillations of the nonhydrostatic pressure.
This test is the 3D analogue of [test2](../test2_temporal_oscillation_modes/instability.c).
The same temporal oscillation modes appear and grow over time with the
default solver. Below the nonhydrostatic pressure is shown in an x-z
slice at the centerline (y ~ 0) of the domain.
Unlike test 2, we do not use temporal relaxation of phi, but instead increase
theta_H slightly. this adds sufficient damping to stabilize the temporal
oscillations.

![Nonhydrostatic pressure fluctuations escalates with theta_H=0.5](instability3D/movie.mp4)

![Nonhydrostatic pressure fluctuations, with theta_H=0.505](instability3D-new/movie.mp4)

*/
//#include <omp.h>
#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"

#include "layered/remap.h"
#include "layered/perfs.h"
#include "test/stokes.h"
#include <sys/stat.h>
#if !NH
scalar w, phi;
#endif
double k_ = 2.*pi, h_ = 0.45, g_ = 1., ak = 0.1;
double RE = 40000.;
#define T0  100.
double slope_angle = 45.; // in degrees
double depth_top = -0.1;
double depth_bottom = -0.30;
double slopex = 0.3; //+- distance from center
double slopedir = 45.;



int main()
{
  origin (-L0/2.,-L0/2.);
  N = 64;
  nl = 25;
  G = g_;
  dt = 0.000001;
  nu = 0.;
  #ifdef implicitness
  theta_H = implicitness;
  #endif
  fprintf(stdout, "Theta_H set to: %g\n", theta_H);
  TOLERANCE = 1e-6 [*];
  DT = 10.;
  run();
}


/* constant flow*/
event init (i = 0)
{
  #if !NH // if NH is not defined:
    assert (nl > 0);
    w = new scalar[nl];
    phi = new scalar[nl];
    reset ({w, phi}, 0.);
  #endif
    CFL_H = 0.25;
    CFL = 0.25;
    
  foreach() {
    double degprmet = 30. [-1];
    double slope = atan((slope_angle-degprmet*sqrt(sq(x)+sq(y)))*3.14159/180.);
    zb[] = -0.5 + x*slope*cos(slopedir*3.14159/180.)+y*slope*sin(slopedir*3.14159/180.);// + (0.005 [1])*noise();
  }
  foreach() {
    double H = 0. - zb[];
    //double z = zb[];
    foreach_layer() {
      h[] = H/nl;
      //z += h[]/2.;
      u.x[] = 0.;
      u.y[] = 0.;
      w[] = 0.;
      //z += h[]/2.;
    }
  }
  boundary(all);
  
}

event logfile (i+=10)
{
  double umax = 0, phimax = 0, phimin = 0;
  foreach (serial) {
    double zc = zb[];
    foreach_layer() {
      double utem = sqrt(sq(w[])+sq(u.x[]));
      if (utem > umax)
        umax = utem;
      if (phi[] > phimax)
        phimax = phi[];
      if (phi[] < phimin)
        phimin = phi[];


    }
  }
  printf ("%g %g %g %g\n", t, umax, phimax, phimin);
}

/**
* x-z slice movie of the nonhydrostatic pressure at the centerline (y ~ 0)
*/
#if 1
event movie (i+=10; t <= T0)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp,
	     "set term pngcairo font ',9' size 800,320\n"
	     "set size ratio -1\n"
	     "set pm3d map interpolate 4,4\n"
	     "unset key\n"
	     "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1,"
	     " 0.25 0 0.5647 1, 0.375 0.05882 1 0.9333,"
	     " 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0,"
	     " 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	     "set cbrange [-1e-6:1e-6]\n"
	     "set cblabel 'phi'\n");
  static int frame = 0;
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "splot [%g:%g][-0.55:0.05] '-' u 1:2:3\n",
	   frame++, t, X0, X0 + L0);
  foreach (serial) {
    if (y > 0 && y < Delta) {
      double z = zb[];
      fprintf (fp, "%g %g %g\n", x, z, phi[]);
      foreach_layer() {
        z += h[];
        fprintf (fp, "%g %g %g\n", x, z, phi[]);
      }
      fprintf (fp, "\n");
    }
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}

event moviemaker (t = end)
{
  system ("rm -f movie.mp4 && "
	  "ffmpeg -y -r 25 -f image2 -pattern_type glob -i 'plot*.png' "
	  "-vcodec libx264 -vf format=yuv420p -movflags +faststart "
	  "movie.mp4 2> /dev/null && "
	  "rm -f plot*.png");
}
#endif

