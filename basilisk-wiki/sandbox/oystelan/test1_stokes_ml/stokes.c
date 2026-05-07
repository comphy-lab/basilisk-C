/**
# Breaking Stokes wave (modified)

This is a modified version of the of the original [stokes.c]() that
tests the stability of the solver for really steep breaking waves.

A very steep, third-order Stokes wave is unstable and breaks, which
makes the nonhydrostatic solver struggle due to a small slope term
inconsistency in the pressure solver.

![Animation of the free-surface, original implementation](stokes/movie.mp4)

By modifying the rhs of the pressure equation, much more stable
simulations can be achieved.

![Animation of the free-surface, patched nh.h](stokes-new/movie.mp4)

The evolution of the energy also shows the instability of the
original version.

~~~gnuplot Evolution of the total energy
set xlabel 'Time'
set ylabel 'Total energy'
plot 'log' u 1:($2+$3) w l t 'Original version', \
     '../stokes-new/log' u 1:($2+$3) w l t 'Patched version'
~~~

The number of iterations of the non-hydrostatic solver is also better (lower).

~~~gnuplot Number of iterations of the non-hydrostatic solver
set ylabel '# iterations'
plot 'perfs' u 1:3 w l t 'Original version', \
     '../stokes-new/perfs' u 1:3 w l t 'Patched version'
~~~

Default solver (0) uses the  original implementation. When switched on
(1), the patched nonhydrostatic pressure solver is used. */

// #define switch_solver 0 // 0 = original implementation, 1 = new impl.


/**
VTK output can optionally be switched on, for easier inspection of
the resulting velocity and pressure field. */

#define vtk_output 0


#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#if switch_solver
# include "../multilayer_stability_patch/nh_patched.h"
#else
# include "layered/nh.h"
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

//#define T0  (k_*L0/sqrt(g_*k_))

#define T0 10.

int main()
{
  origin (-L0/2.);
  periodic (right);
  N = 128;
  nl = 45;  
  G = g_;
  nu = 1./RE;
  TOLERANCE = 1e-5 [*];
  CFL_H = 1;
  theta_H = 0.51;
  #if switch_solver
  theta_R = 0.; // 0.002; // this is not necessary when theta_H > 0.5
  #endif
  run();
}

event init (i = 0)
{
  //  geometric_beta (0.4, true);
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
  fprintf (stderr, "%g %g %g\n", t/(k_/sqrt(g_*k_)), ke/2., g_*gpe + 0.125);
}


#if vtk_output
event output_domain_vtk(t = 0; t += 0.01) {
  fprintf(stdout, "domain vtk output at step: %d, time: %.2f \n", i, t);
  static int j = 0;
  char name[100];
  sprintf(name, "domain_%.6i.vts", j++);
  fprintf(stdout, "written to: %s\n", name);
  FILE* fp = fopen(name, "w");
  output_vts_bin_all_layers(fp,(scalar *){phi}, NULL, false, true);
  fprintf(stdout,"ghosts: %d, bghosts: %d\n", GHOSTS, BGHOSTS);
  fclose(fp);
}
#endif

#if 1
event movie (t += 0.025)
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
	     "set cbrange [-0.5:0.5]\n"
	     "set cblabel 'u.x'\n");
  static int frame = 0;
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "splot [%g:%g][-0.2:0.15] '-' u 1:2:3\n",
	   frame++, t/(k_/sqrt(g_*k_)), X0, X0 + L0);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g\n", x, z, u.x[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g\n", x, z, u.x[]);
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
