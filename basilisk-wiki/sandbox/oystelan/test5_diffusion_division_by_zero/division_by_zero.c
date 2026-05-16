/**
# A small pool of water where the seabed is sloping 45 degrees.
The slope penetrates the surface, so that we have dry cells.

![Initial condition at t=0](division_by_zero/snapshot.png)

This testcase will not run if viscocity nu > 0, due to lack of protection
against division by zero in layered/diffusion.h. 

A modified version of diffusion.h (named diffusion-new.h) includes
extensive division-by-zero protection and runs fine with with the
nonhydrostatic layered solver for all values of nu.


*/

#include "grid/multigrid1D.h"

#if switch_solver
#include "hydro_new.h"
#else
#include "layered/hydro.h"
#endif
#include "layered/nh.h"



#include "layered/remap.h"
#include "layered/perfs.h"

#define T0  100.
double slope_angle = 45.; // in degrees
double depth_top = -0.1 [1];
double depth_bottom = -0.3 [1];
double slopex = 0.3 [1]; //+- distance from center

int main()
{
  origin (-L0/2.);
  periodic (right);
  N = 128;
  nl = 11;
  theta_H = 0.505;
  CFL_H = 1.;
  TOLERANCE = 1e-5 [*];
  dtmax = 0.001;
  nu = 1.E-6;
  DT = 0.1;
  run();
}

/* constant flow*/
event init (i = 0)
{
  double slope = atan(slope_angle*3.14159/180.);
  foreach() {
    zb[] = x < slopex ? (((x-slopex)*slope) > (depth_bottom-depth_top) ? depth_top + 
                         ((x-slopex)*slope) : (x < -slopex ? depth_top : ((x+slopex)*slope) < (depth_top-depth_bottom)
                                               ? (depth_top - ((x+slopex)*slope)) : depth_bottom)) : depth_top;
    // fprintf(stdout,"x: %g, zb: %g\n", x, zb[]);
    double H = max(0.,-0.12 - zb[]);
    double z = zb[];
    foreach_layer() {
      h[] = H/nl;
      z += h[]/2.;
      u.x[] = 0.;
      w[] = 0.;
      z += h[]/2.;
    }
  }
}

#if 1
event logfile (i += 10; t <= T0)
{
  double umax = 0, phimax = 0, phimin = 0;
  foreach (reduction(max:umax) reduction(max:phimax) reduction(min:phimin) ) {
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
  fprintf (stderr, "%g %g %g %g\n", t, umax, phimax, phimin);
  assert (!isnan(umax) && !isnan(phimax) && !isnan(phimin));
}
#endif

event snapshot (i = 0)
{
  FILE * fp = popen ("gnuplot", "w");
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
	   "set cblabel 'phi'\n"
	   "set output 'snapshot.png'\n"
	   "set title 't = %.2f'\n"
	   "splot [%g:%g][-0.32:0.05] '-' u 1:2:3\n",
	   t, X0, X0 + L0);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g\n", x, z, phi[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g\n", x, z, phi[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n");
  pclose (fp);
}




