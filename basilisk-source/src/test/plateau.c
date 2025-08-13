/**
# Scalings for Plateau–Rayleigh pinchoff

An inviscid cylinder of dense liquid is unstable under surface
tension forces. The initial growth of the perturbation is described
by the Rayleigh--Plateau linear stability theory. Close to pinchoff,
non-linear effects cannot be neglected and fully non-linear
solutions must be sought. Using simple similarity arguments, one can
predict that the minimum radius of the deformed cylinder should tend
toward zero as $(t_0-t)^{2/3}$ while the axial velocity should
diverge as $1/(t_0-t)^{1/3}$, with $t_0$ the time of pinchoff.

This test case verifies that these scalings can be recovered using
an axisymmetric VOF calculation.

The figure below illustrates the initial growth, pinchoff and
satellite drop formation.

~~~gnuplot Interfaces at times 0.2, 0.6, breakup and 0.8
set size ratio -1
unset key
unset xtics
unset ytics
unset border
array s[4] = [0, 0.6, 1.4, 2.3]
array q[4] = ["0.2","0.6","0.75585","0.8"]
plot for [j=1:4] for [x0=0:2:2] for [i=-1:1:2] 'prof-'.q[j] u (x0+$1):(i*$2 - s[j]) w l lc -1, \
     for [j=1:4] for [x0=0:2:2] for [i=-1:1:2] 'prof-'.q[j] u (x0+1-$1):(i*$2 - s[j]) w l lc -1
~~~

The animation below shows how adaptivity is used to track the high
curvatures and short timescales close to pinchoff. Up to 18 levels of
refinement are used to capture roughly four orders of magnitude in
characteristic spatial scales.

![Mesh and interface evolution](plateau/movie.mp4)(width="80%")

The scalings for the minimum radius and maximum velocity are given
in the figures below, together with the
theoretical fits. The fit is excellent for at least four orders of
magnitude in timescale. The departures from the power laws close to
pinchoff are due to saturation of the spatial resolution (the
minimum value on the y-axis of the first figure is the grid size
$1/2^{18}$).

~~~gnuplot Evolution of the minimum radius
reset
# We define the breakup time as the time when the axial velocity is maximum
t0="`awk '{ if (NF == 3 && $3 > max) { max = $3; t0 = $1; } } END{ print t0 }' < log`"
set key spacing 1.5
set grid
set xlabel 't_0 - t'
set ylabel 'r_{min}'
set logscale
set format x '10^{%L}'
set format y '10^{%L}'
fit [1e-7:1e-3] a*x**(2./3.) 'log' u (t0 - $1):2 via a
plot [1e-8:][1./2**18:]'log' u (t0 - $1):2 ps 0.5 t '', a*x**(2./3.) t 'x^{2/3}'
~~~

~~~gnuplot Evolution of the maximum axial velocity
set ylabel 'u_{max}'
fit [1e-7:1e-3] a*x**(-1./3.) 'log' u (t0 - $1):3 via a
plot [1e-8:]'log' u (t0 - $1):3 ps 0.5 t '', a*x**(-1./3.) t 'x^{-1/3}'
~~~

For a more detailed description see [Popinet, 2009](/sandbox/popinet/README#popinet2009) and
[Popinet & Antkowiak, 2011](/sandbox/popinet/README#popinet2011d).

## Setup

We solve the axisymmetric, incompressible, variable-density,
Navier--Stokes equations with two phases and surface tension. */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

const int maxlevel = 18;

int main()
{
  origin (-0.5);
  f.sigma = 1.;
  rho1 = 1.;
  rho2 = 1e-2;

  /**
  We consider the inviscid case. Viscosity could be added easily with
  something like: */
  
  // mu1 = ..., mu2 = ...;
  
  run();
}

/**
The initial conditions are a perturbed cylinder with a relatively
large amplitude of deformation (10%). */
  
event init (t = 0) {
  double k = pi, R = 0.2;
  fraction (f, R*(1. + 0.1*sin(k*x)) - y);
}

/**
We allocate a field which will be used to store the position of the
interface relative to the axis of symmetry. */

scalar Y[];

/**
We log the minimum of this position as a function of time as well as
the maximum axial velocity. */

event logfile (i += 5)
{
  position (f, Y, {0,1});
  fprintf (stderr, "%.12f %.12f %.12f\n", t, statsf(Y).min, normf(u.x).max);
}

/**
We generate interface profiles and an animation. */

const double tpinch = 0.75626;

event profiles (t = {0.2, 0.6, tpinch, 0.8})
{
  char name[80];
  sprintf (name, "prof-%g", t);
  FILE * fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);
}

#include "view.h"

event movie (t = 0.6; i += 5; t <= tpinch)
{
  view (fov = 30, near = 0.01, far = 1000,
	tx = -0.111, tz = -0.4,
	width = 1024, height = 680);
  squares (color = "level", min = 6, max = maxlevel, spread = -1);
  draw_vof (c = "f");
  mirror ({0,-1}) {
    squares (color = "level", min = 6, max = maxlevel, spread = -1);
    draw_vof (c = "f");
  }
  save ("movie.mp4");
}

/**
The mesh is adapted "manually" so that the axisymmetric radius of
deformation of the interface is alway described by at least 5 grid
points, down to a maximum level of refinement of maxlevel. */

event adapt (i++)
{

  /**
  We check whether the column is broken. */
  
  position (f, Y, {0,1});
  static bool broken = false;
  if (!broken)
    broken = statsf(Y).min < 1./(1 << maxlevel);

  /**
  The refinement uses maxlevel levels before breakup and 10 after. */

  const double eps = 1e-6;
  refine (level < (broken ? 10 : maxlevel) &&
	  f[] > eps && f[] < 1. - eps &&
	  Delta > Y[]/5.);

  /**
  Cells which do not contain the interface (or which are at a level
  larger than 10 after breakup) are "unrefined". */
  
  unrefine (f[] <= eps || f[] >= 1. - eps ||
	    (broken && level > 10));
}

/**
## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/plateau.html)
* [3D plateau example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/plateau.html)
*/
