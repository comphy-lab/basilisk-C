/**
# VOF tracer around a rising bubble

A bubble rises due to buoyancy. The fluid around the bubble carries a
"VOF tracer" with a constant concentration, which should remain
constant.

This is not the case as illustrated in the figure below.

![Concentration field at $t = 0.5$](vof-tracer/c.png)

The deviations from a constant are quite large.

~~~gnuplot Maximum, minimum and average deviations of the concentration (relative to one)
set xlabel 'Time'
set key top left
plot 'log' u 1:3 w l t 'maximum', '' u 1:4 w l t 'average', \
     '' u 1:2 w l t 'minimum'     
~~~
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/perfs.h"
#include "view.h"

scalar T[];

int main (int argc, char * argv[])
{
  origin (-L0/2., -L0/4.);
  rho2 = 1./100.;
  mu1 = 1e-6;
  mu2 = mu1/rho1*rho2;
  f.sigma = 0.5;
  f.tracers = (scalar *){T};
  TOLERANCE = 1e-4;
  run();
}

event acceleration (i++)
{
  face vector av = a;
  foreach_face(y)
    av.y[] -= 9.8;
}

event init (t = 0)
{
  T.gradient = zero; // minmod2;
  fraction (f, sq(x) + sq(y) - sq(0.125));
  foreach()
    T[] = f[];
  boundary ({T});
}

event logfile (i++)
{
  scalar c[];
  foreach()
    c[] = f[] > 1e-10 ? T[]/f[] : 1.;
  stats s = statsf (c);
  fprintf (stderr, "%g %g %g %g\n", t, s.min - 1., s.max - 1.,
	   s.sum/s.volume - 1.);
}

event end (t = 0.5)
{
  scalar c[];
  foreach()
    c[] = f[] > 1e-10 ? T[]/f[] : 1.;
  view (fov = 23.279, tx = 0.0432061, ty = -0.208829);
  box();
  draw_vof ("f");
  squares ("c", spread = -1);
  save ("c.png");

  // this fails
  stats s = statsf (c);
  assert (fabs(s.max - 1.) < 1e-6 && fabs(s.min - 1.) < 1e-6 &&
	  fabs(s.sum/s.volume - 1.) < 1e-6);
}

#if 0 // for debugging
event dumping (i += 10; t <= 0.5)
{
  scalar c[];
  foreach()
    c[] = f[] > 1e-10 ? T[]/f[] : 1.;
  dump();
}
#endif
