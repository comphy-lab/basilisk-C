/**
# Collapsing cavity makes waves

![Animation of free surface, vertical velocity (left), 
  horizontal velocity (right).](cavity/movie.mp4)(width=602,height=333)

~~~gnuplot Evolution of kinetic, potential and total energies.
plot 'log' u 1:2 w l t 'kinetic', \
     '' u 1:(2.32-$3) w l t 'potential', \
     '' u 1:($2+2.32-$3) w l t 'total'
~~~
*/

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"
#include "view.h"
#include "navier-stokes/perfs.h"

u.n[top] = neumann(0);
p[top] = dirichlet(0);

#define MAXLEVEL 7

int main()
{
  size (5.);
  origin (-L0/8.);
  rho2 = 0.001;
  mu1 = 0.0001;
  mu2 = rho2*mu1/rho1;
  G.x = -1;
  N = 1 << MAXLEVEL;
  run();
}

event init (i = 0) {
  fraction (f, sq(x) + sq(y) - sq(0.5));
  foreach()
    if (x > 0.)
      f[] = 0.;
  boundary ({f});
}

event logfile (i++) {
  double ke = 0., pe = 0.;
  foreach(reduction(+:ke) reduction(+:pe)) {
    double r = rho(f[]);
    ke += r*sq(norm(u))/2.*dv();
    pe += r*G.x*x*dv();
  }
  fprintf (stderr, "%g %g %g\n", t, ke, pe);
}

event outputs (i += 100) {
  dump();
}

event movie (t += 0.03; t <= 20.) {
#if 1  
  view (fov = 21.0715, quat = {0,0,-0.707,0.707},
	tx = -0.000575189, ty = 0.0132879, bg = {1,1,1},
	width = 1208, height = 666, samples = 4);
  box (notics = true);
  draw_vof ("f", filled = -1, fc = {1,1,1});
  squares ("u.x", linear = true);
  mirror (n = {0,1}) {
    draw_vof ("f", filled = -1, fc = {1,1,1});
    squares ("u.y", linear = true);
    box (notics = true);
  }
#else
  view (fov = 23.0909, quat = {0,0,-0.707,0.707},
	tx = 0.477677, ty = 0.0336109, bg = {1,1,1},
	width = 1024, height = 1024, samples = 4);
  box();
  draw_vof ("f", filled = -1, fc = {1,1,1});
  squares ("u.x", linear = true);
#endif
  save ("movie.mp4");
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.01,1e-2,1e-2}, MAXLEVEL);
}
#endif
