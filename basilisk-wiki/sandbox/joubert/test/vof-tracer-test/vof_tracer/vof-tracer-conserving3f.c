/**
# VOF tracer test around a compound drop

This test the "VOF tracer" implementation with three phases but no triple point.
A drop of oil containing water inside fall due to gravity. The fluid around the drop carries a
"VOF tracer" with a constant concentration, which should remain
constant. Furthermore we test here the three phase version of the momentemum conserving function.

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
#include "../../../conserving3f_tracer.h"
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
  TOLERANCE = 1e-4;
  run();
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

event acceleration (i++)
{
  face vector av = a;
  foreach_face(y)
    av.y[] -= 9.8;
}

static scalar * interfaces1 = NULL;

event vof (i++) {
  /**
    We allocate three temporary vector fields to store the three components
    of the momentum and set the boundary conditions and prolongation
    functions. */

  vector q1[], q2[], q3[];
  for (scalar s in {q1,q2,q3})
    foreach_dimension()
      s.v.x.i = -1; // not a vector
  for (int i = 0; i < nboundary; i++)
    foreach_dimension() {
      q1.x.boundary[i] = boundary_q1_x;
      q2.x.boundary[i] = boundary_q2_x;
      q3.x.boundary[i] = boundary_q3_x;
    }
#if TREE
  foreach_dimension() {
    q1.x.prolongation = prolongation_q1_x;
    q2.x.prolongation = prolongation_q2_x;
    q3.x.prolongation = prolongation_q3_x;
  }
#endif

  /**
    We split the total momentum $q$ into its three components $q2$,$q3$ and $q1$
    associated with $f2$, $f3$ and $f1$ respectively. */

  foreach()
    foreach_dimension() {
      double fc1 = clamp(f1[],0,1);
      double fc2 = clamp(f2[],0,1);
      double fc3 = clamp(f3[],0,1);
      q2.x[] = fc2*rho2*u.x[];
      q3.x[] = fc3*rho3*u.x[];
      q1.x[] = fc1*rho1*u.x[];
    }
  boundary ((scalar *){q1,q2,q3});

  /**
    We use the same slope-limiting as for the
    velocity field. */

  foreach_dimension() {
    q1.x.gradient = q2.x.gradient = q3.x.gradient = u.x.gradient;
  }

  /**
    Momentum $C2$ is associated with $1 - f$, so we set the *inverse*
    attribute to *true*. We use the same slope-limiting as for the
    velocity field. 
    We associate the transport of $q1$ and $q2$ with $f$ and transport
    all fields consistently using the VOF scheme. */

  f1.tracers = (scalar *){q1,T};
  f2.tracers = (scalar *){q2};
  f3.tracers = (scalar *){q3};
  vof_advection ({f1,f2,f3},i);

  /**
    We recover the advected velocity field using the total momentum and
    the density */

  foreach()
    foreach_dimension()
    u.x[] = (q1.x[] + q2.x[]+ q3.x[])/rho(f1[],f2[],f3[]);
  boundary ((scalar *){u});

  /**
    We set the list of interfaces to NULL so that the default *vof()*
    event does nothing (otherwise we would transport $f$ twice). */

  interfaces1 = interfaces, interfaces = NULL;
}

/**
  We set the list of interfaces back to its default value. */

event tracer_advection (i++) {
  interfaces = interfaces1;
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

#if 1 // for debugging
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
