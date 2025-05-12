/**
# A Dipole-Wall collision using a 4th-order accurate solver

![Vorticity](dip/omg.mp4)

![Level](dip/lev.mp4)
 */

#define NOSLIP_TOP (1)
#include "nsf4t.h"
int initlevel = 9;
double ue = 1e-4;
scalar * tracers = NULL;

int main() {
  periodic (left);
  L0 = 25;
  X0 = Y0 = -L0/2;
  const scalar muc[] = 1e-4;
  nu = muc;
  N = 1 << initlevel;
  run();
}

event logger (i++) {
  if (i == 0)
    printf ("#t i p.i p.nr p2.i p2.nr cells depth speed divmax\n");
  printf ("%g %d %d %d %d %d %ld %d %g %g\n",
	  t, i, mgp.i, mgp.nrelax, mgp2.i, mgp2.nrelax,
	  grid->tn, grid->maxdepth, perf.speed, mgp2.resa);
}

event init (t = 0) {
  scalar omg[];
  scalar psi[];
  psi[top] = dirichlet_top (0);
  psi[bottom] = dirichlet_bottom (0);
  TOLERANCE = 1e-4;
  double xp = -1.124, yp = 0.423;
  foreach() {
    psi[] = 0;
    double omgl = 0;
    foreach_child()
      omgl -= (x - xp)*exp(-sq(x - xp) - sq(y - yp));
    omg[] = omgl/4.;
  }
  poisson (psi, omg);
  boundary ({psi});
  coord f = {-1., 1.};
  foreach_face()
    u.x[] = f.x*(psi[0,1] + psi[-1,1] - psi[0,-1] - psi[-1,-1])/(4.*Delta);
  boundary ((scalar*){u});
  vorticityf (u, omg);
}

event adapt (i++)
 adapt_flow (ue, 99, 1);

event mov (t += 2) {
  scalar omg[];
  vorticityf (u, omg);
  output_ppm (omg, file = "omg.mp4", n = 600,
	      linear = true, map = cool_warm,
	      min = -0.2, max = 0.2);
  foreach()
    omg[] = level;
  output_ppm (omg, file = "lev.mp4", n = 600);
}

event stop (t = 500);
