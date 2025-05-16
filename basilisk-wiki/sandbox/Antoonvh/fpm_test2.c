/**
# Test for point-based differencing

The x-derivative is estimated on a 32x32 grid using an increasing
number of randomly placed points.

~~~gnuplot example point distribution
set size square
plot 'example' palette
~~~

~~~gnuplot Convergence of x-derivative with p-refinement
set logscale xy
set size square
set cbrange [0.5:3.5]
set xr [10:300]
set xlabel 'sqrt(#part)'
set ylabel 'mean error'
plot 'out' palette t 'Mean order', 100*x**(-3) t '3rd order conv.'
~~~

Well done least squares finite point method!
 */
#include "fpm.h"

Particles parts;

int N_from_part (int np) {
  double snp = sqrt(np);
  int Ni = log(snp)/log(2) - 1;
  return 1<<Ni;
}

int main() {
  L0 = 6;
  X0 = Y0 = -L0/2;
  for (int np = 200; np < 5e4; np *= 1.1) {
    N = N_from_part(np); // not used here for p-refinement test!
    init_grid (32);
    parts = new_particles(np);
    foreach_particle_in(parts) {
      p().x = L0/2*noise();
      p().y = L0/2*noise();
      p().s = exp(-sq(p().x) - sq(p().y));
    }
    ref_outdated = true;
    if (np > 1000 && np <= 1100) {
      FILE * fp = fopen ("example", "w");
      foreach_particle_in(parts)
	fprintf (fp, "%g %g %g\n", x, y, p().s);
      fclose(fp);
    }
    
    double error = 0;
    int tot = 0;
    double order_tot = 0;
    foreach(reduction(+:error) reduction(+:tot)) {
      coord cc = {x,y};
      double coefs[10] = {HUGE, HUGE, HUGE, HUGE, HUGE,
			  HUGE, HUGE, HUGE, HUGE, HUGE};
      int order = least_squares_poly_2D(cc, coefs, parts);
      if (coefs[1] < HUGE) {
	error += fabs(coefs[1] + 2*x*exp(-sq(x) - sq(y)));
	tot++;
	order_tot += order;
      }
    }
    printf ("%g %g %g\n", sqrt(np), error/tot, order_tot/tot);
  }
  free_scalar_data (reference);
  free_p();
}
