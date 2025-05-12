/**
# 3D Sessile drop

This is the 3D equivalent of the [2D](sessile.c) test case.

The volume of a [spherical
cap](https://en.wikipedia.org/wiki/Spherical_cap) of radius $R$ and
(contact) angle $\theta$ is
$$
V = \frac{\pi}{3}R^3(2+\cos\theta)(1-\cos\theta)^2
$$
or equivalently
$$
\frac{R}{R_0} = \left(\frac{1}{4}(2+\cos\theta)(1-\cos\theta)^2\right)^{-1/3}
$$
with $R_0$ the equivalent radius of the droplet
$$
R_0 = \left(\frac{3V}{4\pi}\right)^{1/3}
$$
To test this relation, a drop is initialised as a half-sphere
(i.e. the initial contact angle is 90$^\circ$) and the contact angle
is varied between 30$^\circ$ and 150$^\circ$. The drop oscillates and
eventually relaxes to its equilibrium position. The curvature along
the interface is close to constant.

![Relaxation toward a $30^\circ$ contact angle.](sessile3D/movie.mp4)

Note that shallower angles are [not accessible
yet](/src/contact.h). */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "vof.h"

scalar f[], * interfaces = {f};

#include "tension.h"
#include "view.h"

/**
To set the contact angle, we allocate a [height-function
field](/src/heights.h) and set the contact angle boundary condition on
its tangential components. */

double theta0 = 45;

vector h[];
h.t[back] = contact_angle (theta0*pi/180.);
h.r[back] = contact_angle (theta0*pi/180.);

#define MAXLEVEL 5

int main()
{

  /**
  We use a constant viscosity. */

  mu[] = {.1,.1,.1};
  
  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */

  f.height = h;

  /**
  We set the surface tension coefficient and run for the range of
  contact angles. */
  
  f.sigma = 1.;

  N = 1 << MAXLEVEL;
  for (theta0 = 30; theta0 <= 150; theta0 += 30)
    run();
}

/**
The initial drop is a quarter of a sphere. */

event init (t = 0)
{
  fraction (f, - (sq(x) + sq(y) + sq(z) - sq(0.5)));
}

/**
We log statistics on the maximum velocity, curvature and volume. If
the standard deviation of curvature falls below $10^{-2}$, we assume
that the steady shape is reached and we stop the calculation. */

event logfile (i += 10; t <= 10)
{
  scalar kappa[];
  cstats cs = curvature (f, kappa);
  foreach()
    if (f[] <= 1e-3 || f[] >= 1. - 1e-3)
      kappa[] = nodata;
  stats s = statsf (kappa);
  fprintf (fout, "%g %g %g %g %g %g %g %d %d %d %d\n", t, normf(u.x).max,
	   s.min, s.sum/s.volume, s.max, s.stddev, statsf(f).sum,
	   cs.h, cs.f, cs.a, cs.c);
  fflush (fout);
  if (s.stddev < 1e-2)
    return 1; // stops
}

#if 0
event snapshots (i += 10)
{
  scalar kappa[];
  curvature (f, kappa);
  p.nodump = false;
  dump (buffered = true);
}
#endif

/**
At equilibrium we output the (almost constant) radius, volume, maximum
velocity and time. */

event end (t = end)
{
  scalar kappa[];
  curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 4.*statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g %.5g %.3g %.5g\n",
	   N, theta0, R, s.stddev, V, normf(u.x).max, t);
}

/**
We make a movie of the relaxing interface for $\theta = 30^\circ$. We
use symmetries since only a quarter of the drop is simulated. */

event movie (i += 5; t <= 3)
{
  if (theta0 == 30.) {
    view (fov = 26.6776, quat = {0.474458,0.144142,0.234923,0.836017},
	  tx = -0.0137556, ty = -0.00718937, bg = {1,1,1},
	  width = 758, height = 552);
    draw_vof ("f");
    draw_vof ("f", edges = true);
    cells (lc = {1,0,0});
    mirror (n = {1,0,0}) {
      draw_vof ("f");
      draw_vof ("f", edges = true);
      cells (lc = {1,0,0});
    }
    mirror (n = {0,1,0}) {
      draw_vof ("f");
      draw_vof ("f", edges = true);
      cells (lc = {1,0,0});
      mirror (n = {1,0,0}) {
	draw_vof ("f");
	draw_vof ("f", edges = true);
	cells (lc = {1,0,0});
      }
    }
    save ("movie.mp4");
  }
}

/**
We use refinement based on a smooth version of the volume
fraction. This guarantees constant refinement around the interface,
which seems to be necessary to reach balance (this should be
improved). Note however that this is specific to this test case and
should not generally be used in "production" runs, which should work
fine with the default criterion. */

#if TREE
event adapt (i++) {
#if 1
  scalar f1[];
  foreach()
    f1[] = f[];
  adapt_wavelet ({f1}, (double[]){1e-3}, minlevel = 3, maxlevel = MAXLEVEL);
#else
  adapt_wavelet ({f}, (double[]){1e-4}, minlevel = 3, maxlevel = MAXLEVEL);
#endif
}
#endif

/**
We compare $R/R_0$ to the analytical expression.

~~~gnuplot
reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set arrow from 30,1 to 150,1 nohead dt 2
kappa(theta) = 2.*((2. + cos(theta))*(1. - cos(theta))**2/4.)**(1./3.)
R0(V) = (3.*V/(4.*pi))**(1./3.)
set xtics 30,30,150
plot 2./(kappa(x*pi/180.)) t 'analytical', \
     'log' u 2:(2.*$3/R0($5)) pt 7 t 'numerical'
~~~

## See also

* [2D test](sessile.c) 
*/
