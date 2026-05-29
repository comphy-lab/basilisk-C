/**
# Sessile drop

A sessile drop is a drop of liquid at rest on a solid surface. In the
absence of gravity, the shape of the drop is controlled by surface
tension only. An important parameter is the "contact angle" $\theta$ between
the solid surface and the interface. In the absence of gravity, the
drop is hemispherical and it is easy to show that the relation between
the radius of the drop $R$ and its volume $V$ is (for two-dimensional
drops)
$$
V = R^2 (\theta - \sin\theta\cos\theta)
$$

To test this relation, a drop is initialised as a half-disk (i.e. the
initial contact angle is 90$^\circ$) and the contact angle is varied
between 15$^\circ$ and 165$^\circ$. The drop oscillates and eventually relaxes
to its equilibrium position. This equilibrium is exact to within
machine accuracy. The curvature along the interface is constant.

Note that shallower angles are [not accessible yet](/src/contact.h).

~~~gnuplot Equilibrium shapes for $15^\circ \leq \theta \leq 165^\circ$
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
set xrange [-1:1]
set yrange [0:]
plot 'out' w l, '' u (-$1):2 w l lt 1, 0 lt -1
set term pop
~~~
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#if SMALL_ANGLES
# include "contact-small.h"
#else
# include "contact.h"
#endif
#include "vof.h"
#include "tension.h"

scalar f[], * interfaces = {f};

/**
To set the contact angle, we allocate a [height-function
field](/src/heights.h) and set the contact angle boundary condition on
its tangential component. */

vector h[];
double theta0 = 30;
h.t[bottom] = contact_angle (theta0*pi/180.);
#if SMALL_ANGLES
f[bottom] = contact_fraction (theta0*pi/180.);
h.n[bottom] = contact_normal (theta0*pi/180., f);
#endif

int LEVEL = 4;

int main()
{
  size (1 [1]);

  /**
  We use a constant viscosity. */
  
  mu[] = {.1,.1};

  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */

  f.height = h;


  /**
  We set the surface tension coefficient and run for the range of
  contact angles. */
  
  f.sigma = 1.;

  for (LEVEL = 4; LEVEL <= 6; LEVEL++) {
    init_grid (1 << LEVEL);

    for (theta0 = 15; theta0 <= 45; theta0 += 15)
      run();
    for (theta0 = 60; theta0 <= 180; theta0 += 30)
      run();
  }
}

/**
The initial drop is a quarter of a circle. */

event init (t = 0)
{
  fraction (f, - (sq(x) + sq(y) - sq(0.3)));
}

#if 0
event logfile (i++)
{
  fprintf (fout, "%g %g\n", t, normf(u.x).max);
}

event snapshots (t += 1)
{
  p.nodump = false;
  dump();
}
#endif

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u.x,u.y}, {1e-3,1e-2,1e-2}, LEVEL);
}
#endif

/**
At equilibrium (t = 10 seems sufficient), we output the interface
shape and compute the (constant) curvature. */

double exact (double theta) {
  double vol = pi*sq(0.3)/2.;
  return sqrt (vol/(theta - cos(theta)*sin(theta)));
}

event end (t = 3)
{
  if (LEVEL == 6)
    output_facets (f, stdout);
  
  scalar kappa[];
  curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;
  fprintf (stderr, "level %d %g %.5g %.3g\n", LEVEL, theta0, R/sqrt(V/pi), s.stddev);

  double kexact = 1./exact (theta0*pi/180.);
  double k = s.sum / s.volume;
  fprintf (stderr, "convergence %g %g %g\n",
      theta0, 0.3*pow(2., LEVEL)/L0, (k - kexact)/kexact);
}

/**
We compare $R/R_0$ to the analytical expression, with $R_0=\sqrt{V/\pi}$.

~~~gnuplot Default implementation
reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set arrow from 15,1 to 165,1 nohead dt 2
set xtics 15,15,165
plot 1./sqrt(x/180. - sin(x*pi/180.)*cos(x*pi/180.)/pi) t 'analytical', \
  '<grep "level 4" log' u 3:4 pt 7 t 'LEVEL 4', \
  '<grep "level 5" log' u 3:4 pt 7 t 'LEVEL 5', \
  '<grep "level 6" log' u 3:4 pt 7 t 'LEVEL 6'
~~~

~~~gnuplot Default implementation
reset
set xlabel 'Spatial resolution (grid points per drop radius)'
set ylabel 'Relative error on curvature'
set logscale
set xtics (4.8, 9.6, 19.2)
plot [4:25][0.0005:] \
  '<grep "convergence 15 " log' u 3:(abs($4)) w lp pt 7 t '15 degrees', \
  '<grep "convergence 30 " log' u 3:(abs($4)) w lp pt 7 t '30 degrees', \
  '<grep "convergence 45 " log' u 3:(abs($4)) w lp pt 7 t '45 degrees', \
  '<grep "convergence 60 " log' u 3:(abs($4)) w lp pt 7 t '60 degrees', \
  '<grep "convergence 90 " log' u 3:(abs($4)) w lp pt 7 t '90 degrees', \
  '<grep "convergence 120" log' u 3:(abs($4)) w lp pt 7 t '120 degrees', \
  '<grep "convergence 150" log' u 3:(abs($4)) w lp pt 7 t '150 degrees', \
  '<grep "convergence 180" log' u 3:(abs($4)) w lp pt 7 t '180 degrees', \
  x**-2/3. lt 0 t 'x^{-2}', x**-1/3. lt 0 t 'x^{-1}'
~~~

~~~gnuplot Default implementation
reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set arrow from 15,1 to 165,1 nohead dt 2
set xtics 15,15,165
plot 1./sqrt(x/180. - sin(x*pi/180.)*cos(x*pi/180.)/pi) t 'analytical', \
  '<grep "level 4" ../sessile-gerris-small/log' u 3:4 pt 7 t 'LEVEL 4', \
  '<grep "level 5" ../sessile-gerris-small/log' u 3:4 pt 7 t 'LEVEL 5', \
  '<grep "level 6" ../sessile-gerris-small/log' u 3:4 pt 7 t 'LEVEL 6'
~~~

~~~gnuplot Small angles
reset
set xlabel 'Spatial resolution (grid points per drop radius)'
set ylabel 'Relative error on curvature'
set logscale
set xtics (4.8, 9.6, 19.2)
plot [4:25][0.0005:] \
  '<grep "convergence 15 " ../sessile-gerris-small/log' u 3:(abs($4)) w lp pt 7 t '15 degrees', \
  '<grep "convergence 30 " ../sessile-gerris-small/log' u 3:(abs($4)) w lp pt 7 t '30 degrees', \
  '<grep "convergence 45 " ../sessile-gerris-small/log' u 3:(abs($4)) w lp pt 7 t '45 degrees', \
  '<grep "convergence 60 " ../sessile-gerris-small/log' u 3:(abs($4)) w lp pt 7 t '60 degrees', \
  '<grep "convergence 90 " ../sessile-gerris-small/log' u 3:(abs($4)) w lp pt 7 t '90 degrees', \
  '<grep "convergence 120" ../sessile-gerris-small/log' u 3:(abs($4)) w lp pt 7 t '120 degrees', \
  '<grep "convergence 150" ../sessile-gerris-small/log' u 3:(abs($4)) w lp pt 7 t '150 degrees', \
  '<grep "convergence 180" ../sessile-gerris-small/log' u 3:(abs($4)) w lp pt 7 t '180 degrees', \
  x**-2/3. lt 0 t 'x^{-2}', x**-1/3. lt 0 t 'x^{-1}'
~~~


## See also

* [Similar test with
   Gerris](https://gerris.dalembert.upmc.fr/gerris/tests/tests/sessile.html)
*/
