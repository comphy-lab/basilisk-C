#include "grid/octree.h"
#include "embed.h"
#include "perlin.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
double z0 = 0.05;
#define Z_0 (z0)
#include "rough-embed.h"
#include "navier-stokes/perfs.h"
#include "../profile5c.h"
#include "view.h"
#include "lambda2.h"
double Re = 25000;

double H = 100; // Height of U damp layer with speed U;
coord U = {20., 0 ,1.}; // U.y = 0! 
double muv = 0;
double A0 = 10;

#if DNS
u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
u.r[embed] = dirichlet (0.);
face vector muc[];
#endif

int maxlevel = 9;
double ue = 0.5;

int main(int argc, char ** argv) {
  if (argc > 1)
    maxlevel = atoi(argv[1]);
  if (argc > 2)
    ue = atof(argv[2]);
#if DNS
  if (argc > 3)
    Re = atof(argv[3]);
  muv = U.x*H/Re;
  mu = muc;
#endif
  L0 = 2.*H;
  periodic (left);
  periodic (back);
  run();
}

event init (t = 0) {
  double nx = 6, ny = 1;
  srand (0);
  init_perlin (nx, ny);
  vertex scalar phi[];
  refine (fabs (y - A0) < A0 && level < maxlevel);
  foreach_vertex() 
    phi[] = 5*A0*perlin (x, z + 0.1, nx, ny) + y - A0;
  boundary ({phi});
  fractions (phi, cs, fs);
  output_ppm (cs, file = "cs.png", n = 300);
  view (fov = 30, tx = -0.15, ty = -0.25, phi = 0.3, theta = 0.5);
  draw_vof ("cs", "fs", fc = {sin(pid()),cos(pid()),sin(pid())});
  save ("surface.png");
  foreach() {
    if (y > A0) {
      foreach_dimension() 
      u.x[] = min(U.x*(sqrt((y - A0)/H) + .1), U.x);
    }
  }
  boundary ({cs, fs, u});
}
#if DNS
event properties (i++) {
  foreach_face()
    muc.x[] = muv*fs.x[];
  boundary((scalar*){muc});
}
#endif

event profs (t += 1) {
  char fname[99];
  sprintf (fname, "prof%g", t);
  profile ((scalar*){u}, fname, ym = A0, h = Y0 + L0, n = 30);
}

event forcing (i++) {
  double Hy = A0 + H;
  foreach() {
    if (y > Hy) {
      double damp = ((y - Hy)/(H));
      foreach_dimension() {
	u.x[] += damp*(U.x - u.x[])*dt/5.;
      }
    }
  }
  boundary ((scalar*){u});
}

event adapt (i++) {
  cs.prolongation = refine_injection;
  boundary ({cs});
  adapt_wavelet ({cs, u}, (double[]){0.001, ue, ue, ue}, maxlevel, 1.);
  cs.prolongation = fraction_refine;
  boundary ({cs});
}

event mov (t += 0.5) {
scalar m[];
  scalar l2[];
  foreach()
    m[] = cs[] - 0.5;
  output_ppm (u.x, file = "ux.mp4", mask = m, n = 300., min = -15, max = U.x);
  foreach()
    m[] = level;
  output_ppm (m, file = "lev.mp4", n = 300., min = 1, max = maxlevel);
  view (fov = 30, tx = -0.25, ty = -0.15, phi = 0.3, theta = 0.5);
  draw_vof ("cs", "fs");
  lambda2( u, l2);
  vorticity (u, m);
  isosurface ("l2", -0.1, fc = {0.8, 0.2, 0.8});
  cells(alpha = 0);
  squares ("m", alpha = 0, min = -.5, max = .5, map = cool_warm);
  save ("mov.mp4");
}

event stop (t = 200) {
  dump();
}
