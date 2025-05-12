/**
# VOF tracer test around a compound drop

This test the vOF tracer implementation with three phases but no triple point.
A drop of oil containing water inside fall due to gravity. The fluid around the drop carries a
"VOF tracer" with a constant concentration, which should remain
constant.

This is the case as illustrated in the figure below with the last version of Basilisk.

![Concentration field at $t = 0.3$](vof-tracer3f/c.png)

The deviations from a constant are low and pass the test.

~~~gnuplot maximimum, average and minimum deviation of the concentration value in function of time.
reset
set xlabel 'time [s]'
set ylabel 'concentration-1 quantity []'
plot 'log' u 1:3 w l t 'maximum',\
'log' u 1:4 w l t 'average',\
'log' u 1:2 w l t 'minimum'
~~~
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "../../../three-phase.h"
#include "tension.h"
#include "navier-stokes/perfs.h"
#include "view.h"

scalar T[];

int main (int argc, char * argv[])
{
  origin (-L0/2., 0);
  rho1 = 1/100.;
  rho2 = 1.;
  rho3 = 0.5;
  mu1 = 1e-6;
  mu3 = 2e-6;
  mu2 = mu1/rho1*rho2;
  mu3 = 2.*mu2;
  f1.sigma = (0.5+0.5-0.5)/2.;
  f2.sigma = (0.5-0.5+0.5)/2.;
  f3.sigma = (-0.5+0.5+0.5)/2.;
  f1.tracers = (scalar *){T};
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
  fraction (f1, sq(x) + sq(y-0.7) - sq(0.125));
  fraction (f2, -sq(x) - sq(y-0.7) + sq(0.0625));
  fraction (f3, (-sq(x) - sq(y-0.7) + sq(0.125))*(sq(x) + sq(y-0.7) - sq(0.0625)));
  foreach() {
    T[] = f1[];
  }
  boundary ({T});
}

event logfile (i++)
{
  scalar c[];
  foreach()
    c[] = f1[] > 1e-12 ? T[]/f1[] : 1.;
  stats s = statsf (c);
  fprintf (stderr, "%g %g %g %g\n", t, s.min - 1., s.max - 1.,
	   s.sum/s.volume - 1.);
}

event end (t = 0.3)
{
  scalar c[];
  foreach()
    c[] = f1[] > 1e-12 ? T[]/f1[] : 1.;
  view (fov = 24, tx = 0.0120017, ty = -0.477666);

  box();
  draw_vof ("f1");
  draw_vof ("f2");
  squares ("c", spread = -1);
  save ("c.png");

  // this was failing before the patch
  stats s = statsf (c);
  assert (fabs(s.max - 1.) < 1e-6 && fabs(s.min - 1.) < 1e-6 &&
	  fabs(s.sum/s.volume - 1.) < 1e-6);
}

#if 0 // for debugging
event dumping (t +=0.1; t <= 0.3)
{
  scalar c[];
  foreach()
    c[] = f1[] > 1e-12 ? T[]/f1[] : 1.;
  char name[80];
  sprintf (name, "snap-%g", t);
  dump(file = name);
}
#endif
