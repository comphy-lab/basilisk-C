/**
# Sessile drop on an embedded boundary

This test case is very close to the [sessile
drop](/src/test/sessile.c) test but uses an embedded boundary rather
than the boundary of the domain.

~~~gnuplot Equilibrium shapes for $15^\circ \leq \theta \leq 165^\circ$
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
set xrange [-1.6:1.6]
set yrange [0:]
plot 'out' w l, '' u (-$1):2 w l lt 1, 0 lt -1
set term pop
~~~
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "contact-embed.h"

double theta0;

int main()
{
  size (2.);

  /**
  We shift the bottom boundary. */

  origin (0, - 0.26);
  
  /**
  We use a constant viscosity. */

  mu1 = mu2 = 0.1;
  
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  We vary the contact_angle. */

  for (theta0 = 15; theta0 <= 165; theta0 += 15) {
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
}

event init (t = 0)
{

  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = y;
  boundary ({phi});
  fractions (phi, cs, fs);
  fraction (f, - (sq(x) + sq(y) - sq(0.5)));
}

event logfile (i++; t <= 20)
{

  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-6)
    return true;
}

/**
Given the radius of curvature $R$ and the volume of the droplet $V$,
this function returns the equivalent contact angle $\theta$ verifying
the equilibrium solution
$$
V = R^2(\theta - \sin\theta\cos\theta)
$$
*/

double equivalent_contact_angle (double R, double V)
{
  double x0 = 0., x1 = pi;
  while (x1 - x0 > 1e-4) {
    double x = (x1 + x0)/2.;
    double f = V - sq(R)*(x - sin(x)*cos(x));
    if (f > 0.)
      x0 = x;
    else
      x1 = x;
  }
  return (x0 + x1)/2.;
}

event end (t = end)
{

  /**
  At the end, we output the equilibrium shape. */
  
  output_facets (f, stdout);

  /**
  We compute the curvature only in full cells. */
  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;
  
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g %.4g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi);
}

/**
We compare $R/R_0$ to the analytical expression, with $R_0=\sqrt{V/\pi}$.

The accuracy is not as good (yet) as that of the [sessile
test](/src/test/sessile.c).

~~~gnuplot
reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set arrow from 15,1 to 165,1 nohead dt 2
set xtics 15,15,165
plot 1./sqrt(x/180. - sin(x*pi/180.)*cos(x*pi/180.)/pi) t 'analytical', \
  'log' u 2:3 pt 7 t 'numerical'
~~~

Another way to display the same result is to compare the "apparent
contact angle" with the imposed contact angle.

~~~gnuplot
reset
set xlabel 'Imposed contact angle (degrees)'
set ylabel 'Apparent contact angle (degrees)'
unset key
set xtics 15,15,165
set ytics 15,15,165
plot 'log' u 2:5 pt 7, x
~~~

## See also

* [Sessile drop on an inclined boundary](sessile-inclined.c)
* [Sessile drop on a domain boundary](/src/test/sessile.c)
*/
