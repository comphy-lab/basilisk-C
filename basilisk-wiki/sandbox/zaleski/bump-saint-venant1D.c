/**

# Test code

*/

#include "grid/cartesian1D.h"
#include "saint-venant.h"

int LEVEL = 9;

int main (int argc, char * argv[])
{
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  DT = 1e-1;
  run();
}

event init (i = 0)
{
  foreach()
    {
      h[] = 0.1*L0 + L0*exp(-200.*(x*x + y*y)/(L0*L0));
      u.x[] = 0.;
    }
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (ferr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

event outputfile (t <= 0.4; t += 0.2) {
  foreach()
    {
      // if (h[] < 0.)
      //   return 1; // stops if and when rho becomes unphysical. 
      fprintf (stdout, "%g %g %g\n", x, h[],h[]*u.x[]);
    }
  fprintf (stderr, "\n");
}

/**

# Results
~~~gnuplot  Evolution of the density and momentum. The results are comparable to that obtained when running with the Saint-Venant solver but the MacCormack scheme has oscillations.
plot 'out' u 1:2 w l t "density rho", 'out' u 1:3 w l t "momentum q.x"
# set term @PNG
# set output "bump-saint-venant1D.png"
# replot
~~~

*/
