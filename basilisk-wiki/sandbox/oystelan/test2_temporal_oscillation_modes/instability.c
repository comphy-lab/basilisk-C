/**
# Water at rest over a sloping bottom with relatively steep slopes

This test shows the temporal oscillation modes which over time builds.
Water starts to move and eventually it crashes with the default
solver. 

![Nonhydrostatic pressure fluctuations escalates with the original
 solver](instability/movie.mp4)

If we apply a small degree of relaxation however, the growth
of the temporal mode is stopped and kept at sub-tolerance levels.

![Nonhydrostatic pressure fluctuations, stabilized with
 relaxation](instability-new/movie.mp4)

This also shows on the maximum norm of the velocity.

~~~gnuplot Evolution of the maximum norm of the velocity
set xlabel 'Time'
set ylabel 'Maximum of the norm of velocity'
set logscale y
set yrange [1e-7:]
set key top left
plot 'log' u 1:2 w l t 'Original solver', \
     '../instability-new/log' u 1:2 w l t 'Patched solver'
~~~

**Note that setting theta_H = 0.51 or larger is enough to stabilise
the original version.** */

#include "grid/multigrid1D.h"
#include "layered/hydro.h"

// #define switch_solver 0 // 0 = original implementation, 1 = applied 2% relaxation on phi.
#define output_vtk 0 // write vts files that shows the velocity of each cell

#if switch_solver
# include "../multilayer_stability_patch/nh_patched.h"
#else
# include "layered/nh.h"
#endif

#include "layered/remap.h"
#include "layered/perfs.h"
#if output_vtk
# include "../output_vts_multilayer.h"
# include <sys/stat.h>
#endif

#define T0  200.
double slope_angle = 45.; // in degrees
double depth_top = -0.1;
double depth_bottom = -0.30;
double slopex = 0.3; //+- distance from center

int main()
{
  origin (-L0/2.);
  periodic (right);
  N = 128;
  nl = 11;
  theta_H = 0.500; // setting this to 0.51 stabilises both versions, without theta_R
  CFL_H = 1.;
  TOLERANCE = 1e-5 [*];
  #if switch_solver
  theta_R = 0.02;
  #endif
  dtmax = 0.001;
  DT = 10.;
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
    double H = 0. - zb[];
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
}
#endif

/**
* output vts files which can be viewed in paraview
*/
#if output_vtk
event output_domain_vtk(i+=10; t <= T0) {
    //fprintf(stdout, "domain vtk output at step: %d, time: %.2f \n", i, t);
    struct stat st = { 0 };
    if (stat("vtk", &st) == -1) {
        mkdir("vtk", 0700);
    }
    static int j = 0;
    char name[100];
    sprintf(name, "%s/domain_%.6i.vts", "vtk", j++);
    FILE* fp = fopen(name, "w");
    output_vts_bin_all_layers(fp,(scalar *){phi}, NULL, false, false);
    fclose(fp);
    //fprintf(stdout, "done\n");
}
#endif

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
	   "splot [%g:%g][-0.32:0.05] '-' u 1:2:3\n",
	   frame++, t, X0, X0 + L0);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g\n", x, z, phi[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g\n", x, z, phi[]);
    }
    fprintf (fp, "\n");
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

