/**
# Poisson equation "solver"


~~~gnuplot Source and solution
set xr [-12:12]
set size ratio -1
set cbrange [-1:1]
plot 'out' u ($1-6):2:3 palette pt 5, 'out' u ($1+6):2:4 palette pt 5
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
	p1.x = x + Delta*noise()/3.;
	p1.y = y + Delta*noise()/3.;
	p1.s = 0;//exp(-sq(p1.x) - sq(p1.y));
	p1.b = 4*(-1 + sq(p1.x) + sq(p1.y))*exp(-sq(p1.x) - sq(p1.y));
	add_particle (p1, p);
    }
    }
    assign_particles(p, reference);
    
    poisson(p);
    double errx = 0;
    int np = 0;
    foreach_particle_in(p, reduction(+:errx)) {
      np++;
      double coefs[10] = {0};
      coord X = {x,y,z};
      int order = least_squares_poly_2D (X, coefs, p);
      errx += fabs(coefs[1] + 2*x*exp(-sq(x)-sq(y)));
      if (N == 32)
	printf ("%g %g %g %g\n", x, y, p().b, p().s);
    }
    fprintf (stderr, "%g %g\n", sqrt(np), errx/np);
    free_scalar_data(reference);
  }
  
}
	   
