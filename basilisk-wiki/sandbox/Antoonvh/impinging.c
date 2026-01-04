/**
# Impinging jets

Two opposing jets in a reactor can create a splash. This page focusses
on the visualization. The physical scenario borrows from [the
atomizing jet example](/src/examples/atomization.c), but concerns two
opposing jets.

![`Output_ppm()` is neat for slices, but cant do "3D" shading effects](impinging/f.mp4)

![With Bview, a 3D snapshot can be rendered each iteration with little overhead](impinging/mov.mp4)

![Bwatch can visualize refractive objects](impinging/bw.mp4)
 */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "bwatch.h"
#include "tag.h"

#define radius 1./12.
#define length 0.025
#define Re 400
#define SIGMA 3e-4

double uemax = 0.1;
double femax = 0.01;
int maxlevel = 8;

/** Jet injection with a Poissuille-flow profile prevents early atomization
 */
#define RAD (sqrt(sq(y) + sq(z)))
#define UR (-(RAD - radius)*(RAD + radius)/sq(radius)*min(2., 10.*t))

scalar f0[];
u.n[left]  = dirichlet(f0[]*UR);
u.t[left]  = dirichlet(0);
#if dimension > 2
u.r[left]  = dirichlet(0);
#endif
p[left]    = neumann(0);
f[left]    = f0[];

u.n[right]  = dirichlet(-f0[-1]*UR);
u.t[right]  = dirichlet(0);
#if dimension > 2
u.r[right]  = dirichlet(0);
#endif
p[right]    = neumann(0);
f[right]    = f0[-1];

u.n[top] = neumann(0.);
p[top] = dirichlet (0.);
u.n[bottom] = neumann(0.);
p[bottom] = dirichlet (0.);

int main () {
  init_grid (64);
  origin (0, -1.5, -1.5);
  size (3.);
  rho1 = 1., rho2 = 1./27.84;
  mu1 = 2.*radius/Re*rho1, mu2 = 2.*radius/Re*rho2;  
  f.sigma = SIGMA;
  run();
}

event init (t = 0) {
  refine (1.1*radius > RAD && level < maxlevel);
  fraction (f0, sq(radius) - sq(y) - sq(z));
#if TREE
  f0.refine = f0.prolongation = fraction_refine;
#endif
  restriction ({f0}); // for boundary conditions on levels
  foreach() {
    f[] = f0[]*(x < length) + f0[]*(x > (X0 + L0 - length));
    u.x[] = f[]* sign (X0 + L0/2 - x);
  }
  boundary ({f, u.x});
  // Get an image
  system ("wget https://hips.hearstapps.com/hmg-prod.s3.amazonaws.com/images/screen-shot-2020-09-30-at-1-56-23-pm-1601488673.png?resize=768:*");
  system ("convert screen-shot-2020-09-30-at-1-56-23-pm-1601488673.png?resize=768:* " 
	  "-resize 500x500! reactor.png");
}

event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){femax, uemax, uemax, uemax}, maxlevel, 5);
}
/**
# Movie makers
 */

event mov (t += 0.025) {
  output_ppm (f, file = "f.mp4", n = 500);
}

event bviewer (i++) {
  view (tx = -0.5, phi = 0.3, theta = 0.5);
  draw_vof ("f");
  cells();
  save ("mov.mp4");
}

event droplet_cleanup_project (t += 0.01) {
  // Remove fov crap ...
  remove_droplets (f, minsize = 4);
  // ... and remove atomized features
  remove_droplets (f, minsize = -10);
}

event bwatch_images (t += 0.01) {
  static FILE * fp = popen ("ppm2mp4 bw.mp4", "w");
  watch (fov = 2.5, poi = {1.75, 0.25, 0}, O = {6, 3, 5},
	 nx = 600, ny = 500);
  image ("reactor.png", res = 400, alpha = Z0 + L0/10.);
  sphere (R = 30, mat = {.dull = true});
  equiplane (f, vof = true, mat = {.ind = 1.2});
  volume (f, sc = 0.2, mval = 0.5, col = {0, 95, 104});
  store (fp);
  plain();
}

event stop (t = 2.5) {
  return 1;
}
