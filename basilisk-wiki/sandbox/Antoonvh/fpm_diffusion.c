/**
# Diffusion using the finite point method

![Diffusion of a Gaussian pulse](fpm_diffusion/mov_diffusion.mp4)

~~~gnuplot We obtain (approximately) second-order convergence
set logscale xy
set grid
set xr [20:200]
set xlabel 'sqrt(#parts)'
set ylabel 'Mean error'
plot 'out', 0.5*x**(-2) t '2nd order'
~~~

For stable time integration we need to store the estimated Laplacian
for interpolation of the tendency
*/
#define ADD_PART_MEM double s; double lap;
#include "fpm.h"
#include "run.h"
Particles s;

double kappa = 0.5;

double solution (double x, double y, double t) {
  double t0 = 1;
  return exp((-sq(x)-sq(y))/(4*kappa*(t + t0)))/((t + t0));
}

FILE * gnuplotPipe;

int main() {
  L0 = 14;
  X0 = Y0 = -L0/2.;
  //max_particles = 17;
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf(gnuplotPipe,
	  "set term pngcairo size 700,600\n"
	  "set xr [-7: 7]\n"
	  "set yr [-7: 7]\n"
	  "set size square\n"
	  "set key bottom left\n"
	  "set xlabel 'x'\n"
	  "set ylabel 'y'\n"
	  "set key outside\n"
	  "set cbrange [-0.1:1.1]\n");
  
  for (N = 16; N <= 64; N *= 2) {
    DT = sq(L0/(2*N))/(4*kappa);
    init_grid (N);
    ref_outdated = true;
    run();
  }
  pclose (gnuplotPipe);
  system("rm mov_diffusion.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov_diffusion.mp4");
  system("rm plot*");
}

event init (t = 0) {
  
  s = new_particles (0);
  // initialize with 4 particles per cell, Such that there are up to
  // 6 x 6 particles in a stencil (i.e. more than max_particles)
  foreach() {
    foreach_child() {
      particle p;
      p.x = x + noise()*Delta/4.;
      p.y = y + noise()*Delta/4.;
      p.s = solution(p.x, p.y, t);
      add_particle (p, s);
    }
  }
}
/**
## Diffusion using the interpolated Laplacian
 */

event diffusion (i++, last) {
  
  dt = dtnext(DT);
  // Estimate local 
  foreach_particle_in(s) {
    double coef[10] = {0};
    coord pc = {x,y};
    if (least_squares_poly_2D(pc, coef, s) >= 2) {
      if (coef[4] != HUGE && coef[5] != HUGE)
	p().lap = 2*(coef[4] + coef[5]); 
    } else
      p().lap = 0;
  }
  // Shuffle and store
  foreach_particle_in(s) {
    p().z = p().s;
    p().s = p().lap;
  }
  //interpolate laplacian
  foreach_particle_in(s) {
    double coef[10] = {0};
    coord pc = {x,y};
    least_squares_poly_2D(pc, coef, s);   
    p().lap = coef[0]; //Could their difference be a good refinement indicator?
  }
  // Update solution
  foreach_particle_in(s)
    p().s = p().z + dt*kappa*p().lap;
  // Boundary hack
  foreach_particle_in(s)
    if (fabs(x - (X0 + L0/2)) > L0/2.2 || fabs(y - (Y0 + L0/2)) > L0/2.2)
      p().s = 0;
}

/**
## Movie
 */
int frame = 0;
event movie (t += 0.05) {
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n"
	  "set title 'Diffusion t = %2g'\n", frame++, t);
  fprintf(gnuplotPipe, "plot '-'  palette pt 5 t 'field'\n");
  foreach_particle_in(s)
    fprintf (gnuplotPipe, "%g %g %g\n",x, y, p().s);
  fprintf(gnuplotPipe, "e\n");
}


event stop (t = 1) {
  double error = 0;
  int tot = 0;
  foreach_particle_in(s, reduction(+:error), reduction(+:tot)) {
    error += fabs(p().s - solution(x, y, t));
    tot++;
  }
  printf ("%d %g\n", 2*N, error/tot);
  free_scalar_data (reference);
  return 1;
  
}
