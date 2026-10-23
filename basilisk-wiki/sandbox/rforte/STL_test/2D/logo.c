#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

double L = 0.376;   
double U0 = 50.;
double MU0 = 1e-3;

int maxlevel = 12; 
int minlevel = 5;
double uemax = 0.05;

scalar d[];
face vector muc[];

u.n[left]  = dirichlet (U0);
p[left]    = neumann (0.);
pf[left]   = neumann (0.);

u.n[right] = neumann (0.);
p[right]   = dirichlet (0.);
pf[right]  = dirichlet (0.);

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

int main() {
  size (L);
  origin (-L/2., -L/2.);
  init_grid (1 << minlevel);

  mu = muc;
  run();
}

event init (t = 0) {
  
  restore ("restart");

  solid (cs, fs, (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.);
  fractions_cleanup (cs, fs);

  foreach()
    u.x[] = cs[] > 0. ? U0 : 0.;
}

event properties (i++) {
  foreach_face()
    muc.x[] = MU0*fm.x[];
}

event logfile (i++)
  fprintf (stderr, "%d %g %g %ld\n", i, t, dt, grid->tn);

event movie (t += 2e-4) {
  scalar omega[];
  vorticity (u, omega);
  view (width = 800, height = 800);
  clear();
  draw_vof ("cs", "fs", filled = -1);
  squares ("omega", spread = 20, linear = true);
  save ("vorticity.mp4");
  
  clear();
  view (width = 1640, height = 802);
  cells();
  draw_vof ("cs", "fs", lw = 1.5, lc = {1,0,0});
  box ();
  save ("cells.mp4");
}

event adapt (i++)
  adapt_wavelet ({cs, u}, (double[]){1e-2, uemax, uemax}, maxlevel, minlevel);
solid (cs, fs, (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.);

event end (t = 0.02);


# Video 