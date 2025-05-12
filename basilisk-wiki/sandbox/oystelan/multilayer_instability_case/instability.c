/**
# Modified version of the classic Breaking Stokes wave

Reduced steepness of wave and added a variable bottom topography
where the slope can be adjusted by the parameter slope_angle

instabilities in the no hydrostatic pressure phi occures after
a couple of seconds.

*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "test/stokes.h"
#include "output_vts_multilayer.h"
#include <sys/stat.h>

double k_ = 2.*pi, h_ = 0.45, g_ = 1., ak = 0.1;
double RE = 40000.;
#define T0  (k_*L0/sqrt(g_*k_))
double slope_angle = 15.; // in degrees
double depth_top = -0.45;
double depth_bottom = -0.50;
double slopex = 0.4; //+- distance from center



int main()
{
  origin (-L0/2.);
  periodic (right);
  N = 128;
  nl = 44;
  G = g_;
  nu = 0;
  CFL_H = 0.1;
  TOLERANCE = 1e-6 [*];
  max_slope = 1.; // a bit less dissipative
  DT = 0.001;
  run();
}

event init (i = 0)
{
    double slope = atan(slope_angle*3.14159/180.);
  foreach() {
    fprintf(stdout,"x: %g, ll: %g\n",x, ((x-slopex)*slope));
    zb[] = x < slopex ? (((x-slopex)*slope) > (depth_bottom-depth_top) ? depth_top + ((x-slopex)*slope) : (x < -slopex ? depth_top : ((x+slopex)*slope) < (depth_top-depth_bottom) ? (depth_top - ((x+slopex)*slope)) : depth_bottom)) : depth_top;
    double H = wave(x, 0) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H/nl;
      z += h[]/2.;
      u.x[] = u_x(x, z);
      w[] = u_y(x, z);
      z += h[]/2.;
    }
  }
}

event logfile (i++)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
  printf ("%g %g %g\n", t/(k_/sqrt(g_*k_)), ke/2., g_*gpe + 0.125);
}


/**
* output vts files which can be viewed in paraview
*/
#if 1
event output_domain_vtk(t += T0/4.; t <= 5*T0) {
    fprintf(stdout, "domain vtk output at step: %d, time: %.2f \n", i, t);
    struct stat st = { 0 };
    if (stat("vtk", &st) == -1) {
        mkdir("vtk", 0700);
    }
    static int j = 0;
    char name[100];
    sprintf(name, "%s/domain_%.6i.vts", "vtk", j++);
    fprintf(stdout, "written to: %s\n", name);
    FILE* fp = fopen(name, "w");
    
    output_vts_bin_all_layers_multivar((scalar *){phi},fp);
    fclose(fp);
}
#endif

