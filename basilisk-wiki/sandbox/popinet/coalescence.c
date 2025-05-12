#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "view.h"

scalar f[], * interfaces = {f};

int main()
{
  size (4.);
  origin (-L0/2., -L0/2.);
  const face vector muc[] = {0.01,0.01};
  mu = muc;
  f.sigma = 1.;
  run();
}

event init (t = 0)
{
  fraction (f, max (- (sq(x + 1.) + sq(y) - sq(0.4)),
		    - (sq(x - 1.) + sq(y) - sq(0.5))));
  foreach()
    u.x[] = - sign(x)*f[];
}

event movie (t += 0.04; t <= 6.)
{
  clear();
  squares ("u.x", spread = -1, linear = true);
  draw_vof ("f");
  box();
  save ("movie.mp4");
}
