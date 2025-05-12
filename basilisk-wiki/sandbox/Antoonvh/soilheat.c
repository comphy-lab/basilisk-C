/**
# The temperature in the soil

The interesting and challenging field of atmospheric boundary-layer
 turbulence considers the facinating physics that occur above the
 earth's surface. Underneath, a simple diffusion equation describes
 the evolution of the thermodynamic variable. We use data from the
 excellent study: *missed fog?* of [Izett et
 al. (2019)](https://doi.org/10.1007/s10546-019-00462-3), whom also
 meassured the temperature at the top of the soil (i.e in the
 grass). This is used to force the one-dimensional diffusion problem.

This is the result for the 18-day campaign:

![Profile and grid points](soilheat/mov.mp4)

Even when using adaptive methods, the number of cells is too large to
implement this in global model. If only there was a better way...
 */

#include "grid/bitree.h"
#include "diffusion.h"
#include "run.h"

double * t_top, T_top;        
scalar T[];                   //Temperature field
T[right] = dirichlet (T_top); //boundary condition at the top (i.e. right)
int size_a;

double kappa = 0.000001;      //diffusivity of heat in m^2/s

int main() {
  T.refine = refine_linear;
  X0 = -L0;
  run();
}
/**
## Loading the `data`

We read the data and store the time and the temperature in their
respective arrays. The initial profile is *guessed* to be constant.
*/
event init (t = 0) {
  system ("gunzip -k -f soilheat.temp_cab.bin.gz");
  FILE * fp = fopen ("soilheat.temp_cab.bin", "rb");
  fread (&size_a, 4, 1, fp);
  t_top          = (double*)malloc (size_a*sizeof(double));
  double * times = (double*)malloc (size_a*sizeof(double));
  fread (times, sizeof(double), size_a, fp);
  fread (t_top, sizeof(double), size_a, fp);
  fclose (fp);
  DT = 24*3600*times[size_a - 1]/(double)size_a;  //This is 37sec instead of 30...
  free (times);
  foreach()
    T[] = t_top[0];
}
/**
   Each timestep we update the temperature at the soil top and solve
   for the diffusion problem.
 */

event diff (i++) {
  T_top = t_top[i];
  dt = dtnext (DT);
  const face vector kap[] = {kappa};
  diffusion (T, dt, kap);
}

/**
To capture the evolution proper, a maximum mesh size of $\approx
4\mathrm{mm}$ is used, when its needed. The refinement criterion is
linked to the precision of the meassurements ($\pm 0.1 ^oC$). 
 */

event adapt (i++)
  adapt_wavelet ({T}, (double[]){0.1}, 8);
/**
## Generate a movie

We open a pipe to `gnuplot` and generate a bunch of frames. At the end
of the run we use `ffmpeg` to append these frames into a movie.
*/
FILE * gnuplotPipe;
event init_pipe (t = 0) {
  double minT = HUGE, maxT = -HUGE;
  for (int j = 0; j < size_a ; j++) {
    maxT = max(t_top[j], maxT);
    minT = min(t_top[j], minT);
  }
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf (gnuplotPipe,
           "set term pngcairo\n"
           "set xr [%g: %g]\n"
	   "set yr [-1.0: 0.0]\n"
           "set grid\n"
           "set title 'Soil heat'\n"
	   "set xlabel 'Temperature [^oC]'\n"
	   "set ylabel 'Depth [m]'\n"
	   "set key off\n",
	   minT, maxT);
}

int frame_nr = 0;
event frames (i += 10) {
  fprintf (gnuplotPipe, "set output 'plot%d.png'\n", frame_nr);
  fprintf (gnuplotPipe, "plot '-' w l lw 5, '-' ps 3 pt 2\n");
  foreach_boundary(left)
    fprintf (gnuplotPipe, "%g %g\n",T[], x);
  foreach()
    fprintf (gnuplotPipe, "%g %g\n",T[], x);;
  fprintf (gnuplotPipe, "%g %g\n",T_top, X0 + L0);
  fprintf (gnuplotPipe, "e\n");
  foreach()
    fprintf (gnuplotPipe, "%g %g\n",T[], x);
  fprintf (gnuplotPipe, "e\n");
  frame_nr++;
}

event stop (t = 18*24*3600) {
  pclose (gnuplotPipe);
  system ("rm mov.mp4");
  system ("ffmpeg -loglevel quiet -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system ("rm plot*");
  free (t_top);
  return 1;
}
