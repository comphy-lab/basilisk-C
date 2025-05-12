/**
# VOF tracer around a rising bubble

Adapt mesh test to a tracer attached to 1-f.

A bubble rises due to buoyancy. The fluid around the bubble carries a
"VOF tracer" with a constant concentration, which should remain
constant.

This is the case as illustrated in the figure below with the last version of basilisk.

![Concentration field at $t = 0.5$](vof-tracer_adapt2/c.png)

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

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/perfs.h"
#include "view.h"

scalar T[];
int maxlevel= 7;

int main (int argc, char * argv[])
{
  origin (-L0/2., -L0/4.);
  init_grid (128);
  rho2 = 1./100.;
  mu1 = 1e-6;
  mu2 = mu1/rho1*rho2;
  f.sigma = 0.5;
  T.inverse = true;
  T.refine = T.prolongation = vof_concentration_refine;
  T.c = f;
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
  refine (sq(x) + sq(y) < 1.5*sq(0.125) && level < maxlevel);
  fraction (f, sq(x) + sq(y) - sq(0.125));
  foreach()
    T[] = 1.-f[];
  boundary ({T});
}

event logfile (i++)
{
  scalar c[];
  foreach()
    c[] = (1.-f[]) > 1e-12 ? T[]/(1.-f[]) : 1.;
  stats s = statsf (c);
  fprintf (stderr, "%g %g %g %g\n", t, s.min - 1., s.max - 1.,
	   s.sum/s.volume - 1.);
}

/**
  We adapt the mesh by controlling the error on the volume fraction, 
  velocity field and VOF tracer field. */

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.01,0.01,0.01}, maxlevel);
}

event end (t = 0.5)
{
  scalar c[];
  foreach()
    c[] = (1.-f[]) > 1e-12 ? T[]/(1.-f[]) : 1.;
  boundary ({c});
  view (fov = 23.279, tx = 0.0432061, ty = -0.208829);
  box();
  draw_vof ("f");
  squares ("c", spread = -1);
  save ("c.png");
  // this was failing before the patch
  stats s = statsf (c);
  assert (fabs(s.max - 1.) < 1e-6 && fabs(s.min - 1.) < 1e-6 &&
	  fabs(s.sum/s.volume - 1.) < 1e-6);
}

#if 0 // for debugging
event dumping (t +=0.1; t <= 0.5)
{
  scalar c[];
  foreach()
    c[] = (1.-f[]) > 1e-12 ? T[]/(1.-f[]) : 1.;
  boundary ({c});
  char name[80];
  sprintf (name, "snap-%g", t);
  dump(file = name);
}
#endif
