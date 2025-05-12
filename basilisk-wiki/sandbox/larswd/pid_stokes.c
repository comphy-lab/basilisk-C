/** 
This test is a simple implementation of [pid.h](pid.h) on the [1D breaking stokes wave](http://basilisk.fr/src/test/stokes.c) example of Popinet with PID dampening
*/
 
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "stokes.h"
#include "pid.h"

double k_ = 2.*pi, h_ = 0.5, g_ = 1., ak = 0.35;
double RE = 40000.;
double xmax = 30;
scalar desired_lvl;
#define T0  8*(k_*L0/sqrt(g_*k_))


int main()
{
  origin (-L0/2.);
  periodic (right);
  N = 256;
  nl = 60;
  G = g_;
  nu = 1./RE;
  CFL_H = 1.0;
  max_slope = 1.; // a bit less dissipative
  run();
}

event init (i = 0)
{
  desired_lvl = new scalar[nl];
  double Kp = 1e-2; double Ki = 1e-6; double Kd = 1e-3;
  init_corrector(Kp, Ki, Kd);
  foreach() {
    zb[] = -0.5;
    double H = wave(x, 0) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H/nl;
      z += h[]/2.;
      u.x[] = u_x(x, z);
      w[] = u_y(x, z);
      z += h[]/2.;
      desired_lvl[] = -zb[]/nl;
    }
  }
}

event profiles (t += T0/4.; t <= 2.5*T0) {
  foreach (serial) {
    double H = zb[];
    foreach_layer()
      H += h[];
    fprintf (stderr, "%g %g\n", x, H);
  }
  fprintf (stderr, "\n\n");
}

event pid_control(i++){
  coord x0 = {3.0*L0/4.0, 0, 0};
  coord x1 = {L0, 0,0};
  pid_scalar(desired_lvl, h, x0=x0, x1=x1);
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

event movie (i += 3)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio -1\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][-1.0:0.15]'-' u 1:(-1):2 w filledcu lc 3 t ''",
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
	  "ffmpeg -r 25 -f image2 -i plots/plot%04d.png "
	  "-vcodec libx264 -vf format=yuv420p -movflags +faststart "
	  "movie.mp4 2> /dev/null && "
	  "rm -f plot*.png");
}

event cleanup(t = end){
  delete ((scalar *){desired_lvl});
}