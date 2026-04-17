/**
# Breaking Stokes wave

A steep, third-order Stokes wave is unstable and breaks.

![Animation of the free-surface](stokes/movie.mp4)

The solution obtained using the layered model matches the
Navier-Stokes/VOF solution remarkably well, even after breaking.

~~~gnuplot Wave evolution: layered (left column) and Navier-Stokes/VOF (right column) { width=100% }
unset key
unset xtics
unset ytics
unset border
set multiplot layout 1,2
set size ratio -1
plot for [i = 0:10] 'log' index i u 1:($2-0.15*i) w l lc -1 lt 1
plot for [i = 0:10] '../stokes-ns/log' index i u 1:($2-0.15*i) w l lc -1 lt 1
unset multiplot
~~~

See [Popinet (2020)](/Bibliography#popinet2020) for a more detailed
discussion and [stokes-ns.c]() for the Navier-Stokes/VOF code. 

Modified: to test different solver implementations. new_impl contains changes
in the pressure solver which improves stability substantially.
*/

#define vtk_output 1
#define switch_solver 0 // 0 = original implementation, 1 = new impl.

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#if switch_solver
#include "nh_patched.h"
#else
#include "layered/nh.h"
#endif


#include "layered/remap.h"
#include "layered/perfs.h"
#include "test/stokes.h"
#include "layered/perfs.h"


double k_ = 2.*pi, h_ = 0.5, g_ = 1., ak = 0.45;
double RE = 40000.;
#if vtk_output
#include "../output_vts_multilayer.h"
#endif
#include <sys/types.h>
#include <sys/stat.h>


//#define T0  (k_*L0/sqrt(g_*k_))

#define T0 10.5

int main()
{
  omp_set_num_threads(4);
  origin (-L0/2.);
  periodic (right);
  N = 256;
  nl = 90;  
  G = g_;
  nu = 1./RE;
  TOLERANCE = 1.0e-5 [*];
  theta_H = 0.50;
  #if switch_solver
  theta_R = 0.02;
  #endif
  run();
}

event init (i = 0)
{
  CFL=0.25;
  CFL_H= 0.25;
  //geometric_beta (0.4, true);
  foreach() {
    zb[] = -0.5;
    double H = wave(x, 0) - zb[];
    //double H = - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H/nl;
      z += h[]/2.;
      u.x[] = u_x(x, z);
      w[] = u_y(x, z);
      //u.x[] = 0.1*(z+0.25);
      //w[] = 0.;
      z += h[]/2.;
    }
  }
}


#if 0
event profiles (t += T0/4.; t <= 2.5*T0) {
  foreach (serial) {
    double H = zb[];
    foreach_layer()
      H += h[];
    fprintf (stderr, "%g %g\n", x, H);
  }
  fprintf (stderr, "\n\n");
}
#endif

event logfile (i++; t <= T0)
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


#if vtk_output
event output_domain_vtk(t = 0; t += 0.01) {
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
  output_vts_bin_all_layers(fp,(scalar *){phi}, NULL, false, true);
  fprintf(stdout,"ghosts: %d, bghosts: %d\n", GHOSTS, BGHOSTS);
  fclose(fp);
}
#endif

#if 0
event movie (t += 0.01)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio -1\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][-0.1:0.15]'-' u 1:(-1):2 w filledcu lc 3 t ''",
	   i/3, t/(k_/sqrt(g_*k_)), X0, X0 + L0);
  fprintf (fp, "\n");
  foreach (serial) {
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