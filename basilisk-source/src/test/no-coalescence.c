/**
# Preventing coalescence of VOF tracers

See [no-coalescence.h](/src/no-coalescence.h).

![Animation of the VOF fields: initially only one (red) field is used, then a second (blue) field is introduced to prevent coalescence.](no-coalescence/movie.mp4)

![Animation of the velocity field and interfaces](no-coalescence/velocity.mp4)
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "view.h"
#include "two-phase.h"
#include "no-coalescence.h"

int main()
{
  size (4.);
  origin (-L0/2., -L0/2.);
  mu[] = {0.01,0.01};
  f.sigma = 1.;
  run();
}

event init (t = 0)
{
  fraction (f, max (- (sq(x + 1.) + sq(y) - sq(0.4)),
		    - (sq(x - 1.) + sq(y) - sq(0.5))));
  foreach()
    u.x[] = - 2.*sign(x)*f[];
}

event movie (t += 0.08; t <= 3.)
{
  box();
  int index = 0;
  float * colors[] = {(float[]){1,0,0},(float[]){0,0,1}};
  for (scalar s in interfaces)
    draw_vof (s.name, fc = colors[index++], filled = 1, min = 0, max = 3, lw = 2.);
  save ("movie.mp4");
  clear();

  box();
  for (scalar s in interfaces)
    draw_vof (s.name, fc = {0.13,0.47,0.77}, min = 0, max = 3, lw = 2.);
  squares ("u.x");
  save ("velocity.mp4");
  clear();

  fprintf (stderr, "%g %d\n", t, list_len (interfaces));
}
