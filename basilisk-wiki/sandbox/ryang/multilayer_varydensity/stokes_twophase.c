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

The total energy increases slightly initially (which is not correct)
and decreases linearly after breaking (which is better).

~~~gnuplot Evolution of the energies
reset
set xlabel 'Time'
set key bottom left
plot 'out' u 1:2 w l t 'Kinetic', '' u 1:(($2+$3)/2.) w l t 'Total/2', '' u 1:3 w l t 'Potential', \
     '../stokes-ns/out' u 1:2 w l t 'Kinetic (NS)', '' u 1:(($2+$3)/2.) w l t 'Total/2 (NS)', \
     '' u 1:3 w l t 'Potential (NS)'
~~~
*/

#include "grid/multigrid1D.h"

/**
## Variable density option

We setup a 1/1000 density ratio between the top half and bottom half
(which matches that of the Navier-Stokes/VOF solver). */

#if TWOPHASE
# define NL (nl/2)
# define rho(l) ((l) < NL ? 1.[-3,0,1] : 1./1000.)
# define G_PHI G
#else
# define NL nl
#endif
# define G_ETA 0

/**
## Common setup */

#include "my_hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "test/stokes.h"

double k_ = 2.*pi, h_ = 0.5, g_ = 1., ak = 0.1;
double RE = 40000.;
#define T0  (k_*L0/sqrt(g_*k_))

int writefields_1D (double t, const char *suffix) {
    char filename1[50], filename2[50], filename3[50], filename4[50];
    sprintf (filename1, "field/ux_%s_t%g", suffix, t);
    sprintf (filename2, "field/uz_%s_t%g", suffix, t);
    sprintf (filename3, "field/h_%s_t%g", suffix, t);
    sprintf (filename4, "field/phi_%s_t%g", suffix, t);
    FILE * fux = fopen (filename1, "w");
    FILE * fuz = fopen (filename2, "w");
    FILE * fh = fopen (filename3, "w");
    FILE * fphi = fopen (filename4, "w");
    foreach() {
      foreach_layer() {
        fprintf (fux, "%g %g \n", x, u.x[]);
        fprintf (fuz, "%g %g \n", x, w[]);
        fprintf (fh, "%g %g \n", x, h[]);
        fprintf (fphi, "%g %g \n", x, phi[]);
      }
    }
    fflush (fux);
    fclose (fux);
    fflush (fuz);
    fclose (fuz);
    fflush (fh);
    fclose (fh);
    fflush (fphi);
    fclose (fphi);
    return 0;
}

int main()
{
  origin (-L0/2.);
  periodic (right);
  N = 256;
  nl = 60;
  G = g_;
  nu = 1./RE;
  CFL_H = 0.1;
  max_slope = 1.; // a bit less dissipative
  
  /**
  We setup the variable-density multilayer solver, with a rigid lid. */
  
#if TWOPHASE
  nl *= 2;
  theta_H = 1.;
  rigid = false;
#endif // TWOPHASE
  run();
}

event init (i = 0)
{
  foreach() {
    zb[] = -0.5;
    double H = wave(x, 0) - zb[];
    double z = zb[];
    foreach_layer() {
      if (point.l < NL){
        h[] = H/NL;
        z += h[]/2.;
        u.x[] = u_x(x, z);
        w[] = u_y(x, z);
        z += h[]/2.;
      }
      else {
        h[] = (1. - H)/NL;
        u.x[] = 0.;
        w[] = 0.;
      }
    }
  }
}

event profiles (t += T0/4.; t <= 8.*T0) {
  foreach (serial) {
    double H = zb[];
    foreach_layer()
      if (point.l < NL)
	H += h[];
    fprintf (stderr, "%g %g\n", x, H);
  }
  fprintf (stderr, "\n\n");
}

event logfile (i++)
{
  double kew = 0., gpew = 0., kea = 0., gpea = 0.;
  foreach (reduction(+:kew) reduction(+:gpew) reduction(+:kea) reduction(+:gpea)) {
    double zc = zb[];
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	      norm2 += sq(u.x[]);
      
      if (point.l < NL){
        kew += norm2*h[]*dv();
        gpew += (zc + h[]/2.)*h[]*dv();
        zc += h[];
      }
      else{
        kea += norm2*h[]*dv();
        gpea += (zc + h[]/2.)*h[]*dv();
        zc += h[];  
      }
    }
  }
  static FILE * fp = fopen ("budget.dat", "w");
  fprintf (fp, "%g %g %g %g %g\n", t/(k_/sqrt(g_*k_)), kew/2., g_*gpew + 0.125, kea/2., g_*gpea);
  fflush (fp);
  //printf ("%g %g %g\n", t/(k_/sqrt(g_*k_)), ke/2., g_*gpe + 0.125);
}

event field_log (t += 0.1) {
  char *suffix = "matrix";
  writefields_1D (t, suffix);
}

event movie (i += 3)
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
      if (point.l < NL)
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
