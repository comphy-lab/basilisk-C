/**
# Test for finding nearby particles

Testing the particle-finder with periodic conditions (in left-right
direction), at two levels.

![Finding up to 15 particles](fpm_test/mov.mp4)(loop)

~~~gnuplot We find more particles on the lower levels.
set xlabel 'Level'
set ylabel 'Found particle count'
set xr [2:5]
plot 'out'
~~~
 */
#include "fpm.h"

void plot_nearest (FILE * gnuplotPipe, Particles parts,
		   int * ind, int j, int frame, coord X, int * index_periodic_x) {
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf(gnuplotPipe, "plot '-'  w l lw 2 t 'Nearby', \
          '-' pt 2 t 'All points\n");
  for (int i = 0; i < j; i++) {
    fprintf (gnuplotPipe,"%g %g\n%g %g\n\n",X.x, X.y, pl[parts][ind[i]].x + L0*index_periodic_x[i], pl[parts][ind[i]].y);
  }
  fprintf(gnuplotPipe, "e\n");
  foreach_particle_in(parts)
    fprintf (gnuplotPipe, "%g %g\n",x, y);
  fprintf(gnuplotPipe, "e\n");
}

int main() {
  periodic(left);
  periodic_x = true;
  int index_periodic_x[max_particles];
  L0 = 2;
  X0 = Y0 = -1;
  init_grid(16);
  int np = 50;
  Particles parts = new_particles(np);
  srand(time(NULL));
  foreach_particle_in(parts) {
    p().x = noise();
    p().y = noise();
  }
 
   
  FILE * gnuplotPipe = popen ("gnuplot", "w");
  fprintf(gnuplotPipe,
	  "set term pngcairo\n"
	  "set xr [-1.1: 1.1]\n"
	  "set yr [-1.1: 1.1]\n"
	  "set size square\n"
	  "set key bottom left\n"
	  "set title 'Nearby points (periodic x)'\n"
	  "set xlabel 'x'\n"
	  "set ylabel 'y'\n"
	  "set key outside\n");

  
  
  int frames = 200;
  int frame = 0;
  int levels[2] = {depth(), 3};
  for (int l = 0; l <= 1; l++) {
    int lev = levels[l];
    for (int i = 0; i < frames; i++) {
      fprintf (gnuplotPipe, "set title 'Level = %d'\n", lev);
      coord X = {-0.9*cos(2*pi*frame/frames), 0.9*sin(4*pi*frame/frames) ,0.01};
      int n = 15;
      int ind[n];
      int j = find_nearest_particles(X, n, parts, ind, 2, index_periodic_x, level = lev);
      plot_nearest (gnuplotPipe, parts, ind, j, frame ,X,  index_periodic_x);
      printf ("%g %g\n", lev + noise()/10, j + noise()/5.);
      frame++;
    }
  }
  
  pclose (gnuplotPipe);
  system("rm mov.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system("rm plot*");
}
