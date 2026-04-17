/**
# Water at rest over a sloping bottom with relatively steep slopes.
This test shows the temporal oscillation modes which over time builds,
water starts to move and eventually it crashes with the default
solver. If we apply a small degree of relaxation however, the growth
of the temporal model is stopped and kept at sub-tolerance levels.

*/
//#include <omp.h>
#include "grid/multigrid1D.h"
#include "layered/hydro.h"

#define solvertype 0 // 0 = original implementation, 1 = applied 2% relaxation on phi.
#define output_vtk 1 // write vts files that shows the velocity of each cell

#if solvertype
#include "nh_patched.h"
#else
#include "layered/nh.h"
#endif

#include "layered/remap.h"
#include "layered/perfs.h"
#if output_vtk
#include "../output_vts_multilayer.h"
#endif
#include <sys/stat.h>

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
  theta_H = 0.50;
  TOLERANCE = 1e-5 [*];
  #if solvertype
  theta_R = 0.02;
  #endif
  dtmax = 0.001;
  DT = 10.;
  run();
}

/* constant flow*/
event init (i = 0)
{
    CFL_H = 1.;
    CFL = 0.25;
    double slope = atan(slope_angle*3.14159/180.);
  foreach() {
    
    zb[] = x < slopex ? (((x-slopex)*slope) > (depth_bottom-depth_top) ? depth_top + 
    ((x-slopex)*slope) : (x < -slopex ? depth_top : ((x+slopex)*slope) < (depth_top-depth_bottom)
     ? (depth_top - ((x+slopex)*slope)) : depth_bottom)) : depth_top;
    fprintf(stdout,"x: %g, zb: %g\n",x, zb[]);
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
  boundary(all);
}


#if 1
event logfile (i+=10;t <= T0)
{
  double umax = 0, phimax = 0, phimin = 0;
  foreach () {
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

