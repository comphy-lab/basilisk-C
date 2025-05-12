/**
# Test for Global polynomial remapping

We re-use the wind-driven lake test case as it relies on remapping.

For a dramatic effect, the `rem_frac` is set to 100%.

![Layer heights](wind-driven/mov.mp4)

~~~gnuplot
set yr [-.5:3.5]
set xlabel 'iteration number'
set ylabel 'Remapped cells'
plot 'log' u 1:2 
~~~

~~~gnuplot Some issues remain
reset
unset key
set xlabel 'x'
set ylabel 'z'
scale = 10.
plot [-5:5][0:1]'out' u 1:2:($3*scale):($4*scale) w vectors
~~~

 */

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "gpr.h"

double Re = 10.;
double s = 1./1000.;

double du0;
vector duv[];

int main() {
  L0 = 10.;
  X0 = -L0/2.;
  G = 9.81;
  N = 64;
  nu = sqrt(s*G/Re);
  du0 = sqrt(s*Re*G);
  dut = duv;
  nl = 8;
  run();
}

event init (i = 0) {
  rem_frac = 1;
  foreach() {
    foreach_layer()
      h[] = 1./nl;
    duv.x[] = du0*(1. - pow(2.*x/L0,10));
  }
}

event logger (i++) {
  fprintf (stderr, "%d %d\n", i, nr_remaps);
}

#define uan(z)  (du0*(z)/4.*(3.*(z) - 2.))

event error (t = 1./nu) {
  int i = 0;
  foreach() {
    if (i++ == N/2) {
      double z = zb[], emax = 0.;
      foreach_layer() {
	double e = fabs(u.x[] - uan (z + h[]/2.));
	if (e > emax)
	  emax = e;
	z += h[];
      }
      fprintf (stderr, "# %d %g\n", nl, emax);
    }
  }
}
/**
## movie making
 */
FILE * gnuplotPipe;
event init (t = 0) {
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf(gnuplotPipe,
	  "set term pngcairo\n"
	  "set yr [-5: 5]\n"
	  "set yr [0: 1.1]\n"
	  "set key off\n"
	  "set title 'Lake'\n"
	  "set xlabel 'x'\n"
	  "set ylabel 'h'\n");
}

int frame = 0;
event movie(t += .5){
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf(gnuplotPipe, "plot '-' w l lw 2");
  for (int j = 0; j < nl - 1; j++)
    fprintf(gnuplotPipe, ", '-' w l lw 2");
  fprintf(gnuplotPipe, "\n");
  double dat[grid->tn][nl + 1];
  int n = 0;
  foreach() {
    dat[n][0] = x;
    foreach_layer() {
      z += h[];
      dat[n][point.l + 1] = z;
    }
    n++;
  }
  for (int l = 0; l < nl; l++) {
    for (int j = 0; j < grid->tn; j++) {
      fprintf (gnuplotPipe, "%g %g\n", dat[j][0], dat[j][l + 1]);
      }
    fprintf (gnuplotPipe, "e\n");
  }
  frame++;
}
event stop (t = end){
  system("rm mov.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system("rm plot*");
  return 1;
}

event output (t = end) {
  char name[80];
  sprintf (name, "uprof-%d", nl);
  FILE * fp = fopen (name, "w");
  int i = 0;
  foreach() {
    if (i++ == N/2) {
      double z = zb[];
      foreach_layer()
	fprintf (fp, "%g %g\n", z + h[]/2., u.x[]), z += h[];
    }
    if (nl == 8) {
      double z = zb[];
      foreach_layer()
	printf ("%g %g %g %g\n", x, z + h[]/2., u.x[], w[]), z += h[];
      printf ("\n");
    }
  }
  fclose (fp);
}

