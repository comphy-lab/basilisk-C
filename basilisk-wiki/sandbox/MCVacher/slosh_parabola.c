/**
# Sloshing in a parabolic container (NS)

DNS of the Saint-Venant version of [parabola.c](https://basilisk.fr/src/test/parabola.c). Here the parabola is not as large as in the Saint-Venant case (shallow water hypothesis : $h_0 << \lambda$, not really valid here). However, the linear water profile seems to be still quite valid.

## Simulation
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "navier-stokes/perfs.h"
#include "tag.h"
//#include "harmonic.h"

double h;
double h_1;
double amp;

FILE * fpmax; //

int main() {

  L0=1; //size of the box
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.05;
  mu2 = 0*mu1;
  h=0.31;
  amp=0.8;

  h_1=0.12; //size of the perturbation

  TOLERANCE = 1e-3 [*];

  u.n[embed] = dirichlet(0.);
  u.t[embed] = neumann(0.);

  u.n[top] = neumann(0.);
  p[top] = dirichlet(0.);

  G.y = -9.81;

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
  solid (cs, fs, -(x/amp)*(x/amp)+y-0.02);
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
$(x_G, y_G)$ of the liquid. */

event logfile (i++) {
  double vol = 0., xG = 0., yG = 0.;
  foreach (reduction(+:vol) reduction(+:xG) reduction(+:yG)) {
    double a = f[]*cs[]*dv();
    vol += a, xG += x*a, yG += y*a;
  }
  xG /= vol, yG /= vol;
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%g %g %g\n", t, xG, yG);
  fflush (fpmax);
}

event remove_droplets (i++) {
  remove_droplets (f, threshold=0.05, bubbles=true);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

int isave1 = 1;
event res_save (t += 0.02; t <= 10){

  char name[80];

  sprintf (name, "interface-%d.txt", isave1);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);

  isave1++;
}

event ppm_output (t = 0; t += 0.02; t <= 10) {

  //Note for me : do some tracers only for the water

  scalar m[];
  foreach()
    m[] = cs[] - 0.5;   // hide the solid in the images

  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true, mask = m);

  char name_u_x[80];
  sprintf (name_u_x, "u_x.mp4");
  output_ppm (u.x, file = name_u_x, n = 512, min = -0.3, max = 0.3, linear = true, mask = m);

  char name_u_y[80];
  sprintf (name_u_y, "u_y.mp4");
  output_ppm (u.y, file = name_u_y, n = 512,min = -0.3, max = 0.3, linear = true, mask = m);
}

/**
![Waves](slosh_parabola/f.mp4)
![Vertical velocity](slosh_parabola/u_y.mp4)
![Horizontal velocity](slosh_parabola/u_x.mp4)

~~~gnuplot Horizontal position of the centre of mass of the liquid
set xlabel 't (s)'
set ylabel 'x_G (m)'
set grid
plot 'log.dat' u 1:2 w l lw 2 notitle
~~~
*/