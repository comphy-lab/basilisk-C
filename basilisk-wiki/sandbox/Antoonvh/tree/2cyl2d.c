/**
# Flow around two cylinders in 2D

Following Stephane's setups for the [Von Karman-street
example](/src/examples/karman.c) and the [starting flow around a
cylinder test](/src/tests/starting.c), we also setup a similar
scenario, now using two cylinders.

We use [*automatically scaling* refinement](rc2.c) near the embedded
surface. This is useful when multiple cylinders at various radii
are used.

![Vorticity and tracer particles](2cyl2d/mov2c.mp4)

~~~gnuplot
set xlabel 't'
set ylabel 'C_d'
plot 'log10-10' u 2:3
~~~
 */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "view.h"
#define BVIEW 1
#include "../particles.h"
#include "treegen.h"
#define CYLINDER (sqrt(sq(x) + sq(y - 0.01*R)) - R)

double R = 1, U = 1, Re = 1000, c = 10., C = 2, ue, nu;
int maxlevel = 10;

u.n[left]  = dirichlet (U);
uf.n[left] = U;
u.t[left]  = dirichlet (0);
p[left]    = dirichlet (0.);
pf[left]   = dirichlet (0.);

u.n[right] = neumann (0.);
p[right]   = neumann (0.);
pf[right]  = neumann (0.);

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

FILE * fp;
face vector muc[];
scalar J[], res[];
Branch * branches;


int main() {
  nu = U*R/Re;
  periodic (top);
  L0 = 30;
  X0 = -L0/3.5;
  Y0 = Z0 = -L0/2.;
  mu = muc;
  NITERMIN = 2;
  TOLERANCE = 1e-4;
  char logname[99];
  ue = U/c;
  sprintf (logname, "log%d-%g", maxlevel, c);
  fp = fopen (logname, "w");
  run();
}

void prolongate_ratio (Point point, scalar s) {
  foreach_child() {
    if (s[] != nodata)
      s[] += s[]*Delta;
  }
}

event init (t = 0) {
  res.prolongation = prolongate_ratio;
  nb = 2;
  branches = malloc (nb*sizeof(Branch));
  branches[0].start = (coord){0, 1, -1};
  branches[0].end = (coord){0, 1, 1};
  branches[0].R = R;
  branches[0].parent = 0;
  branches[1].start = (coord){1, -2, -1};
  branches[1].end = (coord){1, -2, 1};
  branches[1].R = 0.35*R;
  branches[1].parent = 1;
  
  foreach()
    u.x[] = U;
  do {
    tree_interface (branches, cs, fs, J);
    foreach() {
      res[] = nodata;
      if (cs[] > 0 && cs[] < 1) {
	res[] = 1/sqrt(branches[(int)(J[] + 0.5)].R/Re);
      }
    }
    boundary ({res, u});
  } while (adapt_wavelet ({res, u}, (double[]){C, ue, ue, ue},
			  maxlevel, 4).nf > 10);
  foreach()
    u.x[] *= (cs[] > 0);
  init_particles_2D_square_grid (10, -5, 0.01, 2);
}
/**
In order to omit issues with inconsitent boundary conditions, we use a
damping layer near the in and out flow boundaries.
 */
event damp (i++) {
  coord Uinf = {U, 0, 0};
  foreach() {
    if (fabs(x - (X0 + L0/2.)) > 4*L0/10.) 
      foreach_dimension()
	u.x[] += dt*(Uinf.x - u.x[])/2.;
  }
  boundary ((scalar*){u});
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*nu; 
  boundary ((scalar*){muc});
}


event adapt (i++) {
  foreach() {
    res[] = nodata;
    if (cs[] > 0 && cs[] < 1) {
      res[] = U/sqrt(branches[(int)(J[] + 0.5)].R*nu);
    }
  }

  adapt_wavelet ((scalar*){res, u}, (double[]){C, ue, ue, ue}, maxlevel, 5);
  unrefine (level > 5 && (x - X0) < L0/10. && (X0 + L0 - x) < L0/10.);
}

event logger (i += 5) {
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  fprintf (fp, "%d %g %g %g %g %g %ld\n",
	   i, t, Fp.x, Fp.y, Fmu.x, Fmu.y, grid->n);
}

event movies (t += 0.4) {
  view (fov = 6, width = 1200, height = 400, tx = -0.25);
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  translate (z = 0.05) {
    draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
    draw_vof ("cs", "fs", lw = 2);
  }
  squares ("omega", min = -2, max = 2, linear = true, map = cool_warm);
  cells();
  char str[99];
  sprintf (str, "Re_R = %g, C = %g, ML = %d", Re, c, maxlevel);
  draw_string (str, 1, lw = 3, lc = {1, 0, 1});
  scatter (loc);
  save ("mov2c.mp4");
}

event stop (t = 150) {
  fclose (fp);
  free (branches);
}
