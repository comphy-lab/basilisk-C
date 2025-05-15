/**
# Test for finding nearby particles

![Finding particles everywhere](fpm_test/mov.mp4)(loop)
 */
#include "fpm.h"

int main() {
  L0 = 2;
  X0 = Y0 = -1;
  init_grid(16);
  Particles parts = new_particles(150);
  foreach_particle_in(parts) {
    p().x = noise();
    p().y = noise();
  }
  FILE * gnuplotPipe = popen ("gnuplot", "w");
  fprintf(gnuplotPipe,
	  "set term pngcairo\n"
	  "set xr [-1: 1]\n"
	  "set yr [-1: 1]\n"
	  "set size square\n"
	  "set key bottom left\n"
	  "set title 'Nearby points'\n"
	  "set xlabel 'x'\n"
	  "set ylabel 'y'\n"
	  "set key outside\n");
  
  int frames = 200;
  for (int frame = 0; frame < frames; frame ++) {
    coord X = {-0.9*cos(2*pi*frame/frames), 0.9*sin(4*pi*frame/frames) ,0.01};
  int n = 15;
  int ind[n];
  int j = find_nearest_particles(X, n, parts, ind, 2);
  
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf(gnuplotPipe, "plot '-'  w l lw 2 t 'Nearby', \
          '-' pt 5 t 'All points\n");
  for (int i = 0; i < j; i++) 
    fprintf (gnuplotPipe,"%g %g\n%g %g\n\n",X.x, X.y, pl[parts][ind[i]].x, pl[parts][ind[i]].y);
  fprintf(gnuplotPipe, "e\n");
  foreach_particle_in(parts)
    fprintf (gnuplotPipe, "%g %g\n",x, y);
  fprintf(gnuplotPipe, "e\n");
  }
  fclose (gnuplotPipe);
  system("rm mov.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system("rm plot*");
}
