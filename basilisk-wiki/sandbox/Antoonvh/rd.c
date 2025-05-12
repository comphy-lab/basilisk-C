/**
# A loosy boundary condition

![](rd/mov.mp4)
 */
#include "grid/multigrid1D.h"
#include "diffusion.h"
#include "run.h"

#define robin(a,b,c) ((dirichlet ((c)*Delta/(2*(b) + (a)*Delta))) \
		      + ((neumann (0))* \
			 ((2*(b) - (a)*Delta)/(2*(b) + (a)*Delta) + 1.)))
#define C_ROBIN (sin(t))

scalar s[];
s[left] = robin (1, 0.1, C_ROBIN);

FILE * gnuplotPipe;

int main() {
  DT = 3e-2;
  TOLERANCE = 1e-5;
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf( gnuplotPipe,
	  "set term pngcairo\n"
	   "set xr [0: 1]\n"
	   "set yr [-1: 1]\n"
	   "set grid\n"
	   "set xlabel 'x'\n"
	   "set ylabel 's'\n");
  run();
}

event differ (i++) {
  dt = dtnext (DT);
  const face vector kappa[] = {.1};
  diffusion (s, dt, kappa);
}

int frame = 0;
event movies (i += 2) {
  fprintf (gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf (gnuplotPipe, "set title 't = %g'\n", t);
  fprintf (gnuplotPipe, "plot '-' w l lw 5 t 's', '' ps 2 t 'c_{robin}'\n");
  foreach()
    fprintf (gnuplotPipe, "%g %g\n",x, s[]);
  fprintf (gnuplotPipe, "e\n");
  fprintf (gnuplotPipe, "0, %g\n e\n", C_ROBIN);
  frame++;
}

event stop (t = 4.*pi) {
  event ("movies");
  pclose (gnuplotPipe);
  system ("rm mov.mp4");
  system ("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system ("rm plot*");
  return 1;
}
