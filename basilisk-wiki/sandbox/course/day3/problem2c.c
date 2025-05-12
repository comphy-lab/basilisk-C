/**
# Influence of the smoothing and resolution on the error of the advection of a function
*/

#include "grid/cartesian1D.h"
#include "utils.h"

double uadv = 1.;

void flux_centered (scalar f, scalar df, double dt)
{
  foreach()
    df[] = (f[1] + f[])/2.*uadv*dt;
}

void flux_upwind (scalar f, scalar df, double dt)
{
  foreach()
    df[] = f[]*uadv*dt;
}

void flux_upwind_second_order (scalar f, scalar df, double dt)
{
  foreach()
    df[] = (3*f[] - f[-1])/2.*uadv*dt;
}

int main() {
  double CFL = 0.1;
  double x0 = 0.3;

  for (int fs = 1; fs <= 16; fs *= 2) {
    for (int n = 5; n <= 5120; n *= 2) {

      double t = 0., tend = 0.2;
      double dt = CFL/n/uadv;
      double sigma = 0.1/fs;

      init_grid (n);

      scalar f[], df[];
      foreach() 
        f[] = erf((x-x0)/sigma);

      //output file
      for (int i = 1; t <= tend; i++) {

        // we obtain the flux (we choose 1st order upwind method)
        //flux_centered (f, df, dt);
        flux_upwind (f, df, dt);
        //flux_upwind_second_order (f, df, dt);

        // numerical solution at t+dt
        foreach()
          f[] -= (df[]-df[-1])/Delta;

        t += dt;
      }

      scalar e[];
      foreach()
        e[] = fabs(f[] - erf((x - (x0 + uadv*t))/sigma));

      printf ("%d %g %g \n", n, sigma, normf(e).max);
    }
    printf(" \n");
  }
}

/**
~~~gnuplot Maximum error as a function of the smoothing factor and resolution
set log xy
set xlabel '1/dx'
set ylabel 'f-f_{exact}'
set cblabel '{/Symbol s}'
p "out" u 1:3:2 not w lp palette pt 7, x**(-1) t 'x^{-1}' w l
~~~ 
 
~~~gnuplot Maximum error as a function of the smoothing factor and resolution (rescaled)
set log xy
set xlabel '{/Symbol s}^2/dx'
set ylabel 'f-f_{exact}'
set cblabel '{/Symbol s}'
p "out" u ($1*$2**2):3:2 not w lp palette pt 7, 0.01*x**(-1) t '0.01 x^{-1}' w l
~~~ 
 */ 