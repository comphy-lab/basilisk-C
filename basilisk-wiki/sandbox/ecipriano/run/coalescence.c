/**
# Droplet Coalescence with Phase Change

We modify the [coalescence.c](/sandbox/popinet/coalescence.c) file using the clsvof approach and including phase change with a simple isothermal evaporation model.

![Evolution of the velocity field](coalescence/movie.mp4)

![Evolution of the species mass fractions](coalescence/species.mp4)
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "two-phase-clsvof.h"
#include "integral.h"
#include "evaporation.h"
#include "species-gradient.h"
#include "view.h"

double Dmix1, Dmix2, YIntVal;

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
uext.n[top] = neumann (0.);
uext.t[top] = neumann (0.);
pext[top] = dirichlet (0.);

u.n[bottom] = neumann (0.);
u.t[bottom] = neumann (0.);
p[bottom] = dirichlet (0.);
uext.n[bottom] = neumann (0.);
uext.t[bottom] = neumann (0.);
pext[bottom] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);
uext.n[left] = neumann (0.);
uext.t[left] = neumann (0.);
pext[left] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);

int main()
{
  size (4. [*]);
  origin (-L0/2. [*], -L0/2. [*]);

  rho1 = 1. [*], rho2 = 1. [*];
  mu1 = 1e-2 [*], mu2 = 1e-2 [*];
  Dmix1 = 0., Dmix2 = 1.e-2;
  YIntVal = 0.6;

  const scalar sigmav[] = 1.;
  d.sigmaf = sigmav;

  run();
}

event init (t = 0)
{
  fraction (f, max (- (sq(x + 1.) + sq(y) - sq(0.5)),
		    - (sq(x - 1.) + sq(y) - sq(0.5))));
  foreach()
    d[] = max (- (sq(x + 1.) + sq(y) - sq(0.5)),
        - (sq(x - 1.) + sq(y) - sq(0.5)));

  foreach() {
    u.x[] = - sign(x)*f[];
    uext.x[] = u.x[];
  }

  foreach() {
    YL[] = f[];
    YG[] = 0.;
    Y[] = YL[] + YG[];
  }
}

event movie (t += 0.04; t <= 6.)
{
  clear();
  squares ("u.x", spread = -1, linear = true);
  draw_vof ("f");
  box();
  save ("movie.mp4");

  clear();
  squares ("Y", min = 0., max = 1., linear = true);
  draw_vof ("f");
  box();
  save ("species.mp4");
}

