/**
# Tidally-induced internal lee waves

This is a modified version of [this example](/src/examples/lee.c),
where the initial and inflow stratification uses an *terrain and
domain following* height. See,

![Evolution of the temperature field](lee2/movie.mp4)

The comments only concern the modifications
*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"

#define drho(T) (1e-3*(T - 13.25)/(8. - 13.25))
/**
The initial condition is refined here */
#define bottom_height(x) (x < 0. ?				\
			  -50. + 35./(1. + pow(x/500.,4)) :	\
			  -100. + 85./(1. + pow(x/500.,4)))
#define zp(z,x) (-100*(z/bottom_height(x))) 
#define T0(z,x) (8. + (13.25 - 8.)*(zp(z,x) + 100.)/100.)

#include "layered/dr.h"
#include "layered/remap.h"
#include "layered/perfs.h"

double nu_H = 0.1;

int main() {
  L0 = 21500;
  X0 = - 6500;
  G = 9.81;
  /** The resolution is reduced for faster computations
   */
  nl = 50;  // 100
  N = 1024; // 2048
  cell_lim = mono_limit;
  DT = 100.;
  nu = 1e-3;
  system ("rm -f plot*.png");
  run();
}

#define M2 (12.*3600. + 25.2*60.)

double Tleft (Point point) {
  double H = 0.;
  double zc = zb[];
  for (int l = - point.l; l < nl - point.l; l++) {
    H += h[0,0,l];
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;
  return T0(zc, x);
}

event init (i = 0) {
  u.n[left]  = dirichlet (0.3*sin(2.*pi*(t/M2)));
  u.n[right] = neumann(0.);
  h[right] = dirichlet(100./nl);
  T[left] = Tleft(point);
  T[right] = Tleft(point);
  /**
     The sill geometry, initial layer depths and initial temperature
     profile. */
  foreach() {
    zb[] = bottom_height(x); 
    double z = zb[];
    foreach_layer() {      
      h[] = - zb[]/nl;
      z += h[]/2.;
      T[] = T0(z,x);
      z += h[]/2.;
      printf ("%g %g\n", z, T[]);
    }
  }
}

event viscous_term (i++) {
  if (nu_H > 0.) {
    scalar d2u[];
    foreach_layer() {
      foreach()
	d2u[] = (u.x[1] + u.x[-1] - 2.*u.x[])/sq(Delta);
      foreach()
	u.x[] += dt*nu_H*d2u[];
#if NH
      foreach()
	d2u[] = (w[1] + w[-1] - 2.*w[])/sq(Delta);
      foreach()
	w[] += dt*nu_H*d2u[];
      boundary ({w});
#endif // NH
    }
    boundary ((scalar *){u});
  }
}

void setup (FILE * fp) {
  fprintf (fp,
#if ISOPYCNAL
	   "set pm3d map corners2color c2\n"
#else
	   "set pm3d map interpolate 2,2\n"
#endif
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set cbrange [8:13.5]\n"
	   "set xlabel 'x (m)'\n"
	   "set ylabel 'depth (m)'\n"
	   "set xrange [-1500:2000]\n"
	   "set yrange [-100:1]\n"
	   );
}

void plot (FILE * fp) {
  fprintf (fp,
	   "set title 't = %.2f M2/8'\n"
	   "sp '-' u 1:2:4\n",
	   t/(M2/8.));
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);  
}

int frame = 0;
event gnuplot (t += M2/512) {
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp, "set term x11 size 1024,300\n");
  if (i == 0)
    setup (fp);
  plot (fp);
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1024,300\n"
	   "set output 'plot%d.png'\n"
           "replot\n", frame++);
}

event figures (t <= M2/2.; t += M2/8.) {
  FILE * fp = popen ("gnuplot 2> /dev/null", "w");  
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1024,300\n"
	   "set output 'T-%g.png'\n", t/(M2/8.));
  setup (fp);
  plot (fp);
}

event moviemaker (t = end) {
  system("rm movie.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y movie.mp4");
}