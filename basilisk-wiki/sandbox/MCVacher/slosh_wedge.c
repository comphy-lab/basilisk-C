/**
# Sloshing in a 90° wedge

We reproduce the sloshing case in a 90° wedge, which resolution is in Lamb's book (*Hydrodynamics*, §258).

## Simulation
*/
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "navier-stokes/perfs.h"
#include "tag.h"

double h;
double h_1;
double H, omega; // liquid depth above the apex, Lamb's frequency

FILE * fpmax;

int main() {

  L0=1;
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.05;
  mu2 = 0*mu1;
  h=0.5;

  h_1=0.12;

  TOLERANCE = 1e-3 [*];

  u.n[embed] = dirichlet(0.);
  u.t[embed] = neumann(0.);

  u.n[top] = neumann(0.);
  p[top] = dirichlet(0.);

  G.y = -9.81;

  H = h - 0.2;            // the apex of the wedge is at y = 0.2
  omega = sqrt (9.81/H);  // see the Theory section

  N=256;
  origin (-L0/2, 0);
  init_grid(N);

  fpmax =  fopen("log.dat", "w");

  run();
}

/**
With `embed`, `f` is the volume fraction of the fluid part of a cell:
it must not be intersected with the solid (this gives $f < 1$ in the cut
cells and spurious bubbles along the wall). It is initialised from the
free surface only. */

event init (t = 0) {
  fraction (f, h_1*x+h-y);
  solid (cs, fs, -fabs(x)+y-0.2);
}

/**
The solid cells are never advected: after each step, the solid cells
next to the fluid take the (cut-cell weighted) average of their fluid
neighbours. This acts roughly as a neutral (about 90°) contact
condition. */

event solid_fraction (i++) {
  foreach()
    if (cs[] == 0.) {
      double fsum = 0., wsum = 0.;
      foreach_neighbor (1)
	if (cs[] > 0.)
	  fsum += cs[]*f[], wsum += cs[];
      if (wsum > 0.)
	f[] = fsum/wsum;
    }
}

/**
The oscillations are monitored through the centre of mass
$(x_G, y_G)$ of the liquid, compared with Lamb's solution for $x_G$
(see below). */

event logfile (i++) {
  double vol = 0., xG = 0., yG = 0.;
  foreach (reduction(+:vol) reduction(+:xG) reduction(+:yG)) {
    double a = f[]*cs[]*dv();
    vol += a, xG += x*a, yG += y*a;
  }
  xG /= vol, yG /= vol;
  double xG_th = 2./3.*h_1*H*cos (omega*t);
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%g %g %g %g\n", t, xG, yG, xG_th);
  fflush (fpmax);
}

event remove_droplets (i++) {
  remove_droplets (f, threshold=0.05, bubbles=true);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

int isave1 = 1;
event res_save (t += 0.02; t <= 25){

  char name[80];

  sprintf (name, "interface-%d.txt", isave1);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);

  isave1++;
}

event ppm_output (t = 0; t += 0.02; t <= 25) {

  scalar m[];
  foreach()
    m[] = cs[] - 0.5;   // hide the solid in the images

  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true, mask = m);

  char name_u_x[80];
  sprintf (name_u_x, "u_x.mp4");
  output_ppm (u.x, file = name_u_x, n = 512, min = -0.8, max = 0.8, linear = true, mask = m);

  char name_u_y[80];
  sprintf (name_u_y, "u_y.mp4");
  output_ppm (u.y, file = name_u_y, n = 512,min = -1, max = 1, linear = true, mask = m);
}

/**
![Waves](slosh_wedge/f.mp4)
![Vertical velocity](slosh_wedge/u_y.mp4)
![Horizontal velocity](slosh_wedge/u_x.mp4)

## Theory (Lamb, §258)

Take the origin at the apex of the wedge, $X = x$, $Y = y - 0.2$, so
that the walls are $Y = \pm X$ and the free surface at rest is
$Y = H$, with $H = h - 0.2$. For an inviscid, irrotational flow, the
potential
$$
\phi = C\, X Y \sin \omega t
$$
is harmonic and satisfies the no-penetration condition on both walls
($\nabla\phi \cdot \mathbf{n} \propto Y - X = 0$ on $Y = X$, and
likewise on $Y = -X$). The linearised free-surface condition
$\partial_t^2\phi + g\,\partial_Y\phi = 0$ on $Y = H$ gives
$$
\omega^2 = \frac{g}{H}.
$$
The free surface stays **plane** and simply rocks about its centre:
$\eta(x,t) = s(t)\, x$ with $s(t) = h_1 \cos \omega t$, which is
exactly our initial condition (released from rest).

To first order in $h_1$, the liquid area is $H^2$ and its first moment
is $\int_{-H}^{H} x\, s x \,\mathrm{d}x = 2 s H^3/3$, so
$$
x_G(t) = \frac{2}{3}\, h_1 H \cos \omega t .
$$
Here $H = 0.3$, so $T = 2\pi\sqrt{H/g} \approx 1.10$ s and the
amplitude of $x_G$ is $0.024$. The theory is inviscid and linear
(the small parameter is the surface slope, here $h_1 = 0.12$), so a slow drift in phase and amplitude is
expected over 25 s.

~~~gnuplot Horizontal position of the centre of mass of the liquid
set xlabel 't (s)'
set ylabel 'x_G (m)'
set grid
plot 'log.dat' u 1:2 w l lw 2 t 'Basilisk', \
     'log.dat' u 1:4 w l dt 2 lw 2 t 'Lamb (linear, inviscid)'
~~~
*/