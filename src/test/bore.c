/**
# Undular bores for the Green-Naghdi equations

This test case was proposed by [Le Métayer et al,
2010](/src/references.bib#lemetayer2010) (section 6.2). This is a dam
break problem described by the (dispersive) [Green-Naghdi
equations](/src/green-naghdi.h) (rather than the non-dispersive
Saint-Venant equations). */

#include "grid/bitree.h"
#include "green-naghdi.h"

/**
The domain is 600 metres long, centered on the origin. The
acceleration of gravity is set to 10 m/s^2^. The problem is
solved in one dimension with 2048 grid points. */

int main()
{
  X0 = -300.;
  L0 = 600.;
  G = 10.;
  N = 2048;
  run();
}

/**
The initial conditions are zero velocity and a jump in fluid depth at
the origin (i.e. dam break conditions). */

event init (i = 0)
{
  foreach() {
    h[] = x < 0. ? 1.8 : 1.;
    u.x[] = 0.;
  }
}

event output (t = 48) {
  foreach()
    fprintf (stdout, "%g %g %g\n", x, h[], u.x[]);
  fprintf (stdout, "\n");
}

/**
At $t = 48$ seconds, the depth and velocity profiles are given
below. They are compared with the numerical solution of the same
problem obtained with the Saint-Venant solver ([bore1.c]()).

The solution consists of localized undular bores superposed onto the
Saint-Venant solution. This demonstrates the robustness of the
numerical scheme.

~~~gnuplot Fluid depth profile at $t = 48$ seconds.
set xlabel 'x'
set ylabel 'z'
set key top left
plot '../bore1/out' w l t 'Saint-Venant', 'out' w l t 'Green-Naghdi'
~~~

~~~gnuplot Velocity profile at $t = 48$ seconds.
set xlabel 'x'
set ylabel 'u'
plot '../bore1/out' u 1:3 w l t 'Saint-Venant', 'out' u 1:3 w l t 'Green-Naghdi'
~~~
*/
