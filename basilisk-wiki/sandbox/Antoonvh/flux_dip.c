/**
# Dynamic flux partitioning example

We consider a net surface flux that is being partitioned between a
solid and a fluid in the following scenario:

![A dipole vortex is targeted at the surface at some moment, after
 the surface has already been receiving some flux. The dipole-induced
 flow creates a very thin thermal boundary which increases the flux
 towards the fluid](flux_dip/mov.mp4)

Next we plot the ratio of the soil flux ($G$) and fluid flux ($H$) 
 
 ~~~gnuplot flux partioning
 set yr [0.08:0.2]
 set xlabel 'time'
 set ylabel 'H/G'
 set grid
 plot 'out' u 1:4 w l lw 2 t 'Flux ratio'
 ~~~
 */

#define FLUXFUN (1)
#include "diffusion-pair.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "view.h"

u.n[embed] = dirichlet (0);
u.t[embed] = dirichlet (0);

scalar s1[], s2[];
scalar * tracers = {s1};

double yv = 1, sv = 0.3;
double yh = -1;

int maxlevel = 8;

face vector muc[];
double muv = 1e-4;

int main() {
  L0 = 4;
  X0 = Y0 = -L0/2;
  mu = muc;
  N = 128;
  DT = .5;
  run();
}

event properties (i++) {
  foreach_face()
    muc.x[] = muv*fs.x[];
}

event init (t = 0) {
  vertex scalar phi[];
  foreach_vertex()
    phi[] = y - yh - sq(x)/10.;
  fractions (phi, cs, fs);
  restriction ({cs, fs});
}

event prepare (t = 49.9) {
  DT = 0.001;
}
event set_dtmax (t = 50) {
  refine (sq(x) + sq(y - yv) < sq(2*sv) && level < maxlevel); 
  scalar omg[], psi[];
  foreach() {
    psi[] = 0;
    omg[] = 20*x*exp(-sq(x/sv)-sq((y - yv)/sv));
  }
  poisson (psi, omg);
  foreach() {
    u.x[] = -cs[]* (psi[0,1] - psi[0,-1])/(2*Delta);
    u.y[] = cs[]*(psi[1] - psi[-1])/(2*Delta);
  }
  DT = 0.1;
}

event tracer_diffusion (i++) {
  kappav1 = 0.001;
  kappav2 = 0.01;
  double rCp1 = 50;
  double rCp2 = 500;
  scalar l1[], l2[];
  foreach() {
    l1[] = rCp1 * cs[];
    l2[] = rCp2 * (1 - cs[]);
  }
  const face vector k1[] = {kappav1, kappav1};
  const face vector k2[] = {kappav2, kappav2};
  diffusion_pair (s1, s2, dt, k1, k2, theta1 = l1, theta2 = l2);
}

event mov (i += 5) {
  scalar omg[];
  vorticity (u, omg);
  view (width = 800, height = 400, fov = 20);
  translate (x = -L0/2) {
    squares ("omg", min = -2, max = 2);
    cells();
  }
  foreach()
    omg[]  = cs[]*s1[] + (1 - cs[])*s2[];
  translate (x = L0/2) {
    squares ("omg", map = cool_warm);
    char str[99];
    sprintf (str, "t = %d", (int)(t + 0.5));
    draw_string (str, size = 30, pos = 2);
    draw_vof ("cs", "fs");
  }
  save ("mov.mp4");
}

event flux_part (i += 5) {
  double H = 0, G = 0;
  foreach (reduction(+:H) reduction(+:G), noauto) {
    if (cs[] > 0 && cs[] < 1.) {
      double V[2];
      embed_flux_pair (point, {s1, s2}, V, false);
      H += V[0]*sq(Delta);
      G += V[1]*sq(Delta);
    }
  }
  printf ("%g %g %g %g\n", t, H, G, H/G);
}

event adapt (i++) {
  adapt_wavelet ({cs, s1, s2, u}, (double[]){1e-4, .1, .1, 1e-2, 1e-2},
		 maxlevel, maxlevel - 3);
}

event stop (t = 200);

