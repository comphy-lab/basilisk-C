
/**
# A bimodal system

We model the evolution of a pdf that statisfies the Fokker-Planck
equation with minimalistic ingredients for a diffusion fueled
bimodality.

The initial Gaussian pdf is placed at the right-hand-side mode.

![Evolution of the pdf in the plotted Drift and diffusion landscape](fpe/mov.mp4)
 */

#include "grid/cartesian1D.h"
#include "fpe.h"

// No drift flux:
rho[left] = 0.; 
rho[right] = 0.;

double D10 = 0;

int main() {
  L0 = 10.;
  X0 = -L0/2;
  DT = 0.001;
  N = 128;
  run();
  /**
  We also displace the "stable point" a bit
  */
  D10 = 0.1;
  run();
}

event init (t = 0) {
  foreach() {
    D2[] = exp(-sq(x)) + 0.01;
    rho[] = exp(-sq(x - 1.));
  }
  foreach_face()
    D1.x[] = -0.1*x + D10;
  stats m = statsf(rho);
  foreach()
    rho[] /= m.sum;
}

scalar rhon[]; 

event stop (t += 1) {
  if (change (rho, rhon) < 1e-4) {
    stats m = statsf (rho);
    printf ("# %g %g\n", t, m.sum);
    if (D10 > 0)
      event ("finalize_movie");
    return 1;
  }
}

event too_long (t = 200){
  if (D10 > 0)
    event ("finalize_movie");
  return 1;
}

/**
## Movie maker code

`Gnuplot` and `FFmpeg` are used to create a movie from the solution
data.
 */

FILE * gnuplotPipe = NULL;

event init (t = 0) {
  if (gnuplotPipe == NULL) {
    gnuplotPipe = popen ("gnuplot", "w");
    fprintf(gnuplotPipe,
            "set term pngcairo\n"
            "set xr [-5: 5]\n"
            "set yr [-0.5: 1.1]\n"
            "set key box top left\n"
            "set grid\n"
            "set title 'PDF evolution'\n"
            "set xlabel 'x'\n"
            "set ylabel 'pdf, D_1, D_2'\n");
  }
}

int frame = 0;
event movie(t += 0.2){
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf(gnuplotPipe, "plot \
          '-' w l lw 5 t 'pdf',			\
          '-' w l lw 2 lc 'black' t 'Drift', \
          '-' w l lw 2 lc 'black' lt ':' t 'Diffusion'\n");
  for (scalar s in {rho, D1.x, D2}) {
    foreach()
      fprintf(gnuplotPipe, "%g %g\n",x, s[]);
    fprintf(gnuplotPipe, "e\n");
  }
  frame++;
}

event finalize_movie(0) {
  system("rm mov.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system("rm plot*");
  pclose (gnuplotPipe);
}

