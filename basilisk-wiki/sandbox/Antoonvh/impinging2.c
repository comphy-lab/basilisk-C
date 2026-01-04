/**
# Hommage to Xiadong Chen et al. (2013)

A [movie of Xiadong Chen](https://www.youtube.com/watch?v=4HEbIRAcwuA) was the inspiration for bview.

Here we marginally reproduce their imfuential results.

![$Re = 200$, $We = 100$, $\rho_1\rho_2^{-1} = 100$, 11 levels of refinement](https://www.antoonvanhooft.nl/media/xiadong_chen.mp4)

![Impinging jets (10 levels)](impinging2/xiadong_chen.mp4)
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
#include "bwatch.h"
#include "tag.h"

#define RAD (sqrt(sq(fabs(y) - dist/2.) + sq(z)))
#define UR (-(RAD - radius)*(RAD + radius)/sq(radius)*min(2., 1e2*t))

// Physical paramters
double dist = 10;          // Distance between jets
double Re = 200, We = 100; // Reynolds and Weber numbers
double angle = pi*30/180.; // Half the angle of impingment

double uj = 1; // Normalized injections veloicty
double radius = 1./2.; // Normalized jet diameter

//Grid paramters
double uemax = 0.05;
double femax = 0.005;
int maxlevel = 10;

// Boundary conditions
// Inlet
scalar f0[];
u.n[left]  = dirichlet(f0[]*uj*cos(angle)*UR);
u.t[left]  = dirichlet(-sign(y)*f0[]*uj*sin(angle)*UR);
#if dimension > 2
u.r[left]  = dirichlet(0);
#endif
p[left]    = neumann(0);
f[left]    = f0[];

u.n[top] = neumann(0.);
p[top] = dirichlet (0.);
u.n[bottom] = neumann(0.);
p[bottom] = dirichlet (0.);

int main () {
  init_grid (64);
  origin (0, -25, -25);
  size (50.);
  rho1 = 100., rho2 = 1.; // Light Water and heavy air (to help numerics)
  // Noting here that with the `conserving` option, the  viscous term is most computationally expensive 
  mu1 = 2.*radius*uj/Re*rho1;
  mu2 = mu1/200.; // (?) Not in paper! From movie ...
  f.sigma = 2*radius*rho1*sq(uj)/We;
  run();
}

event init (t = 0) {
  TOLERANCE = 1e-4;
  refine (1.1*radius > RAD && level < maxlevel && x < 2*radius);
  fraction (f0, radius - RAD);
#if TREE
  f0.refine = f0.prolongation = fraction_refine;
#endif
  restriction ({f0}); // for boundary conditions on levels
  CFL = 0.3;
  boundary ({f, u.x});
}

event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){femax, uemax, uemax, uemax}, maxlevel, 5);
}
/**
# Movie makers
 */

event mov (i += 5) {
  output_ppm (f, file = "f.mp4", n = 500);
}

event bviewer (i += 2) {
  view (tx = -0.5, phi = 0.3, theta = 0.5);
  draw_vof ("f");
  cells();
  save ("mov.mp4");
}

event droplet_cleanup_project (i += 10) {
  // Remove fov crap ...
  remove_droplets (f, minsize = 4);
  // ... and remove atomized features
  remove_droplets (f, minsize = -8);
  foreach() {
    if (x > X0 + 8*L0/9.)
      f[] = 0;
  }
}
#if dimension==3


event bwatch_images (t += 0.05) {
  static FILE * fp = popen ("ppm2mp4 xiadong_chen.mp4", "w");
  //FILE * fp1 = fopen ("test.ppm", "w");
  // Transition from side -> top -> top-side view
  double cam_angle = (1./2.*pi*(tanh((t - 8.)/2.) + 1.)/2. -
		      1./4.*pi*(tanh((t - 20.)/2.) + 1.)/2.);

  watch (fov = 2*dist, poi = {L0/3, 0, 0},
	 O = {L0/3., sin(cam_angle)*L0, cos(cam_angle)*L0},
	 nx = 600, ny = 1000, up = {-1, 1e-9, 1e-9});
  sphere (R = 2*L0, mat = {.dull = true, .n1 = {-1, 0, 0}, .col = {255,255,255}, .col2 = {1,1,1}});
  equiplane (f, vof = true, mat = {.ind = 1.2});
  volume (f, sc = 0.5, mval = 0.5, col = {1, 1, 1});
  store (fp);
  //store (fp1);
  //fclose (fp1);
  plain();
}
#endif
event stop (t = 30.) {
  dump("impinging_angle");
}
/**
## Reference

Chen, Xiaodong, et al. "High-fidelity simulations of impinging jet atomization." Atomization and sprays 23.12 (2013).

*/
