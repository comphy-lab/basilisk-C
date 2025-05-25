/**
# Advection using the finite point method

![Advection of a Gaussian pulse](fpm_diffusion/mov_advection.mp4)
*/
#include "fpm.h"
#include "run.h"
Particles s;

FILE * gnuplotPipe;

coord U = {1, 1};

double solution (double x, double y, double t) {
  return exp(-sq(x-U.x*t + 1) - (sq(y)));
}

int main() {
  periodic (left);
  periodic_x = true;
  periodic (bottom);
  periodic_y = true;
  DT = 0.01;
  init_grid (32);
  L0 = 10;
  X0 = Y0 = -L0/2.;
  run();
}
 
event init (t = 0) {
  s = new_particles(0);
  
  foreach() {
    foreach_child() {
      particle p;
      p.x = x + noise()*Delta/4;
      p.y = y + noise()*Delta/4;
      p.s = solution(p.x, p.y, t);
      add_particle (p, s);;
    }
  }
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf(gnuplotPipe,
	  "set term pngcairo size 700,600\n"
	  "set xr [-5: 5]\n"
	  "set yr [-5: 5]\n"
	  "set size square\n"
	  "set key bottom left\n"
	  "set xlabel 'x'\n"
	  "set ylabel 'y'\n"
	  "set key outside\n"
	  "set cbrange [-0.1:1.1]\n");
}

event advection (i++, last) {
  dt = dtnext(DT);
  foreach_particle_in(s) {
    double coef[10] = {0};
    coord pc = {x,y};
    if (least_squares_poly_2D(pc, coef, s) >= 2) {
      if (coef[1] != HUGE)
	p().z = coef[1]*U.x + coef[2]*U.y; //d/dx
    } else
      p().z = 0;
  }
  foreach_particle_in(s) {
    p().s -= dt*(p().z);
  }
}

int frame = 0;
event movie (t += 0.05) {
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n"
	  "set title 'Advection t = %2g'\n", frame++, t);
  fprintf(gnuplotPipe, "plot '-'  palette pt 5 t 'field'\n");
  foreach_particle_in(s)
    fprintf (gnuplotPipe, "%g %g %g\n",x, y, p().s);
  fprintf(gnuplotPipe, "e\n");
}

event stop (t = 10) {
  smooth_2D(s);
  event("movie");
  pclose (gnuplotPipe);
  system("rm mov.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov_advection.mp4");
  system("rm plot*");
  free_scalar_data(reference);
  return 1;
}
