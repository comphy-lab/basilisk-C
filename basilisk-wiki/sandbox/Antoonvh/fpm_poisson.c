/**
# Poisson equation "solver"


~~~gnuplot Source and solution
set xr [-12:12]
set size ratio -1
set cbrange [-1:1]
plot 'data' u ($1-6):2:3 palette pt 5, 'data' u ($1+6):2:4 palette pt 5
~~~

~~~gnuplot Convergence of the x-derivative
reset
set size square
set logscale x 2
set logscale y
set xr [16:256]
set yr [0.0001:0.01]
set xlabel 'sqrt(#parts)'
set ylabel 'error'
set grid
plot 'log' t  'data', 5*x**-2 t 'second order'
~~~


 */
#include "fpm_poisson.h"

int main() {
  periodic(left);
  periodic(bottom);
  periodic_x = true;
  periodic_y = true;
  L0 = 10;
  X0 = Y0 = -L0/2;
  for (N = 16; N <= 64; N *= 2) {
    init_grid (N);
    Particles p = new_particles (0);
    foreach() {
      foreach_child() {
	particle p1;
#if FPM_BOUNDARY
	p1.bound = 0;
#endif
	p1.x = x + Delta*noise()/3.;
	p1.y = y + Delta*noise()/3.;
	p1.s = exp(-sq(p1.x) - sq(p1.y));
	//p1.b = 4*(-1 + sq(p1.x) + sq(p1.y))*exp(-sq(p1.x) - sq(p1.y));
#if FPM_BOUNDARY
	if (sq(p1.x - 2) + sq(p1.y) < 9)
#endif
	  add_particle (p1, p);
      }
    }
    
    
    #if FPM_BOUNDARY
    int nb = sq(N/5);
    for (int i = 0; i < nb; i++) {
      particle p1;
      p1.bound = 1;
      double theta = 2*pi*i/(double)nb;
      p1.x = 3*cos(theta) + 2;
      p1.y = 3*sin(theta);
      p1.s = exp(-sq(p1.x) - sq(p1.y));
      p1.ds = 0;
      add_particle (p1, p);
    }
    for (int i = 0; i < nb; i++) {
      particle p1;
      p1.bound = 1;
      double theta = 2*pi*(i + 0.5)/(double)nb;
      p1.x = (3+0.5*(L0/N))*cos(theta) + 2;
      p1.y = (3+0.5*(L0/N))*sin(theta);
      p1.s = exp(-sq(p1.x) - sq(p1.y));
      p1.ds = 0;
      add_particle (p1, p);
    }
    for (int i = 0; i < nb; i++) {
      particle p1;
      p1.bound = 1;
      double theta = 2*pi*(i )/(double)nb;
      p1.x = (3+(L0/N))*cos(theta) + 2;
      p1.y = (3+(L0/N))*sin(theta);
      p1.s = exp(-sq(p1.x) - sq(p1.y));
      p1.ds = 0;
      add_particle (p1, p);
    }
#endif
    assign_particles(p, reference);
    foreach_particle_in (p) {
      coord X = {x,y,z};
      double coef[10];
      least_squares_poly_2D (X, coef, p, true);
      p().b = 2*(coef[4]+coef[5]);
    }
    foreach_particle_in (p) {
      p().s = 0;
    }
    poisson(p);
    double errx = 0;
    int np = 0;
    foreach_particle_in(p, reduction(+:errx)) {
      np++;
      double coefs[10] = {0};
      coord X = {x,y,z};
      int order = least_squares_poly_2D (X, coefs, p);
      errx += fabs(coefs[1] + 2*x*exp(-sq(x)-sq(y)));
    }
    
    if (N == 32) {
      FILE * fp = fopen ("data", "w");
      foreach_particle_in(p) 
	fprintf (fp, "%g %g %g %g\n", x, y, p().b, p().s);
      fclose (fp);
      }
    if (N == 64) {
       FILE * fp = fopen ("data_res", "w");
       double max_res = 0;
      foreach_particle_in(p, reduction(max:max_res)) {
#if FPM_BOUNDARY
	// Skipping boundary points 
	if (p().bound == 0) {
#endif
	  coord X = {x, y, z};
	  double coef[10] = {0};
	  least_squares_poly_2D(X, coef, p, true);
	  double lap = 2*(coef[4] + coef[5]);
	  p().res =  p().b - lap;
	  if (fabs(p().res) > max_res)
	    max_res = fabs(p().res);
	  fprintf (fp, "%g %g %g\n", x, y, p().res);
#if FPM_BOUNDARY
	}
#endif
      }
      fclose (fp);
    }
    fprintf (stderr, "%g %g\n", sqrt(np), errx/np);
    free_scalar_data(reference);
  }
  free_p();
  for (scalar s in all)
    printf ("s = %s\n", s.name);
}
	   
