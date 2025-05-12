/**
# Sessile drop on an embedded boundary inclined at 45 degrees

This is a less trivial version of [this test case](sessile.c).

~~~gnuplot Equilibrium shapes for $15^\circ \leq \theta \leq 165^\circ$
set term push
set term @SVG size 480,480
set size ratio -1
unset key
unset xtics
unset ytics
unset border
plot 'out' w l, x + 0.97 with filledcurves x1 lc 'grey'
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
  origin (-1, 0);
  
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
  We define the inclined wall and the initial (half)-circular
  interface. */
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (y - x - 0.97);
  boundary ({phi});
  fractions (phi, cs, fs);
  fraction (f, - (sq(x - 0) + sq(y - 1.) - sq(0.25)));
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
  double R = s.volume/s.sum, V = statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g %.4g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi);
}

/**
We compare $R/R_0$ to the analytical expression, with $R_0=\sqrt{V/\pi}$.

The results are not very good, especially for small contact angles.

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
contact angle" with the imposed contact angle. It is clear that
contact angles smaller than $\approx 45$ degrees or larger than
$\approx 135$ degrees cannot be imposed.

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

* [Sessile drop on a horizontal embedded boundary](sessile.c)
*/
