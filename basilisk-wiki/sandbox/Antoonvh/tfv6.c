/**
# 6th order Finite volume advection

![](tfv6/mov.mp4)
 */

#define RKORDER (4)
#define BGHOSTS 2
#include "grid/bitree.h"
#include "lsrk.h"
#include "higher-order.h"

void advfv6 (scalar * sl, scalar * dsdtl) {
  scalar sa = sl[0];
  scalar sp[];
  reduce_average (sa, sp);
  boundary ({sp});
  foreach() {
    for (scalar dsdt in dsdtl)
      dsdt[] = -(sp[1] - sp[0])/(Delta);
  }
  boundary (dsdtl);
}

scalar sa[];

int main() {
  periodic (left);
  L0 = 20;
  X0 = -L0/2.;
  DT = L0/(1.65*N); 
  run();
}

event init (t = 0) {
  foreach() {
    double a = 0;
    foreach_child() foreach_child()
      a += exp(-sq(x));
    sa[] = a/4.;
  }
}

event advance (i++, last) {
  dt = dtnext (DT);
  A_Time_Step ({sa}, dt, advfv6);
  boundary ({sa});
}

  
#if 1
FILE * gnuplotPipe;
event init_pipe (t = 0) {
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf (gnuplotPipe,
           "set term pngcairo\n"
           "set xr [%g: %g]\n"
	   "set yr [-0.1: 1.1]\n"
           "set grid\n"
	   "set size square\n"
           "set title 'Advection'\n"
	   "set xlabel 'x'\n"
	   "set ylabel 's'\n",
	   X0, X0 + L0);
}

int frame_nr = 0;
event frames (i++) {
  fprintf (gnuplotPipe, "set output 'plot%d.png'\n", frame_nr);
  fprintf (gnuplotPipe, "plot '-' t 'data',"
	   "exp(-(x - %g)**2) lw 2 t 'Analytical'\n", fmod(t + 10, 20) - 10);
  foreach() 
    fprintf (gnuplotPipe, "%g %g\n",x, sa[]);;
  fprintf (gnuplotPipe, "e\n");
  frame_nr++;
}

event stop (t = 100) {
  pclose (gnuplotPipe);
  system ("rm mov.mp4");
  system ("ffmpeg -loglevel quiet -r 25 -f image2 -i plot%d.png \\"
	  "-c:v libx264 -vf format=yuv420p -y mov.mp4");
  system ("rm plot*");
  return 1;
}
#endif
