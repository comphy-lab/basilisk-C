/**
# Collapsing cavity makes waves

![Animation of free surface, vertical velocity (left), 
  horizontal velocity (right).](blow_unstable/movie.mp4)(width=602,height=333)

*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"
#include "curvature.h" // for position()
#include "view.h"
//#include "navier-stokes/perfs.h"

u.n[top] = neumann(0);
p[top] = dirichlet(0);

u.n[right] = dirichlet (- 60.*(y < 0.25));

#define MAXLEVEL 7

int main()
{
  size (5.);
  origin (-L0/2.);
  rho2 = 0.001;
  mu1 = 0.01;
  mu2 = rho2*mu1/rho1;
  G.x = -1;
  N = 1 << MAXLEVEL;
  run();
}

event init (i = 0) {
  foreach()
    f[] = x < 0;
  boundary ({f});
}

/*event logfile (i++) {
  double ke = 0., pe = 0.;
  foreach(reduction(+:ke) reduction(+:pe)) {
    double r = rho(f[]);
    ke += r*sq(norm(u))/2.*dv();
    pe += r*G.x*x*dv();
  }
  scalar X[];
  position (f, X, (coord){1.,0.});
  foreach (reduction(+:c)){
    if (fabs(y - 0.55) <= Delta/2){
      c+= Delta*f[];
    }
  fprintf (stderr, "%g %g %g %g\n", t, ke, pe, c);
}
*/

event outputs (i += 100) {
  dump (buffered = true);
}

event movie (t += 0.03; t <= 5.) {
  view (fov = 21.0715, quat = {0,0,-0.707,0.707},
	tx = -0.000575189, ty = 0.0132879, bg = {1,1,1},
	width = 1208, height = 666, samples = 4);
  box (notics = true);
  draw_vof ("f", filled = 0, fc = {1,1,1});
  squares ("u.x", linear = true);
  mirror (n = {0,1}) {
    draw_vof ("f", filled = -1, fc = {1,1,1});
    squares ("u.y", linear = true);
    box (notics = true);
  }
  save ("movie.mp4");
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.01,1e-2,1e-2}, MAXLEVEL);
}
#endif
