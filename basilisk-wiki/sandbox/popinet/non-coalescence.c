/**
# Non-coalescence of colliding droplets

Two droplets are launched toward one another. This example shows how
using different VOF tracers for each droplet prevents numerical
coalescence. */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "view.h"

scalar f1[], f2[], * interfaces = {f1,f2};

int main()
{
  size (4.);
  origin (-L0/2., -L0/2.);
  const face vector muc[] = {0.01,0.01};
  mu = muc;
  f1.sigma = f2.sigma = 1.;
  run();
}

event init (t = 0)
{
  fraction (f1, - (sq(x + 1.) + sq(y) - sq(0.4)));
  fraction (f2, - (sq(x - 1.) + sq(y) - sq(0.5)));
  foreach()
    u.x[] = f1[] - f2[];
}

event movie (t += 0.04; t <= 6.)
{
  clear();
  squares ("u.x", spread = -1, linear = true);
  draw_vof ("f1");
  draw_vof ("f2");
  box();
  save ("movie.mp4");
}

/**
![Non-coalescence using two VOF tracers. The color field is the
 horizontal component of the velocity
 field.](non-coalescence/movie.mp4)

This can be compared with the [result](coalescence/movie.mp4) for 
a single VOF tracer. The corresponding code is [here](coalescence.c). */
