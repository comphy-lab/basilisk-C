/**
# Buoyant plume

![Animation of free surface, vertical velocity (left), 
  horizontal velocity (right).](plume/movie.mp4)(width=604,height=333)

~~~gnuplot Maximum height as a function of time
plot 'log' u 1:4 w l t ''
~~~
*/

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"

scalar s[];
double rho3 = 0.5;
# define rho(f) (clamp(f,0.,1.)*(clamp(s[],0.,1.)*(rho3 - rho1) + \
				 rho1 - rho2) + rho2)
#include "two-phase.h"

#include "tracer.h"
#include "diffusion.h"
scalar * tracers = {s};

#include "curvature.h" // for position()

#include "view.h"
#include "navier-stokes/perfs.h"

#define MAXLEVEL 7

double radius = 0.1, U0 = 2., tinj = 10.;

u.n[left] = dirichlet (U0*(y < radius && t < tinj));
s[left] = dirichlet(y < radius);

u.n[right] = u.n[] > 0. ? neumann(0) : dirichlet(0);
p[right] = dirichlet(0);

int main()
{
  size (5.);
  origin (-L0/2.);
  rho2 = 0.001;
  mu1 = 0.001;
  mu2 = rho2*mu1/rho1;
  N = 1 << MAXLEVEL;
  const face vector G[] = {-1.,0.};
  a = G;
  run();
}

event tracer_diffusion (i++) {
  const face vector D[] = {0.01,0.01};
  diffusion (s, dt, D);
}

event init (i = 0) {
  foreach()
    f[] = x < 0.;
  boundary ({f});
}

event logfile (i++) {
  double ke = 0., pe = 0.;
  foreach(reduction(+:ke) reduction(+:pe)) {
    double r = rho(f[]);
    ke += r*sq(norm(u))/2.*dv();
    pe -= r*x*dv();
  }
  scalar X[];
  position (f, X, (coord){1.,0.});
  fprintf (stderr, "%g %g %g %g\n", t, ke, pe, statsf(X).max);
}

event outputs (i += 100) {
  dump();
}

event movie (t += 0.03; t <= 20.) {
  view (fov = 21.0715, quat = {0,0,-0.707,0.707},
	tx = -0.000575189, ty = 0.0132879, bg = {1,1,1},
	width = 1208, height = 666, samples = 4);
  box (notics = true);
  draw_vof ("f", filled = -1, fc = {1,1,1});
  squares ("u.x", linear = true);
  mirror (n = {0,1}) {
    draw_vof ("f", filled = -1, fc = {1,1,1});
    squares ("s", linear = true);
    box (notics = true);
  }
  save ("movie.mp4");
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.01,1e-2,1e-2}, MAXLEVEL);
}
#endif

