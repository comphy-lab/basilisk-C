/**
# Von Karman vortex street

![Vorticty](karman/o.mp4)
![grid](karman/l.mp4)
*/
#include "nsrk.h"

u.n[left] = dirichlet(1);
u.n[right] = neumann (0);
p[right] = dirichlet (0);
p2[right] = dirichlet (0);

double R = 0.5, x0 = 15;

int maxlevel = 10;

int main() {
  L0 = 100;
  Y0 = -L0/2;
  N = 128;
  nu = 1./160.;
  run();
}

event init (t = 0) {
  CFL = 0.8;
  foreach_dimension()
    u.x.refine = refine_linear;
  foreach()
    u.x[] = 1;
  boundary ((scalar*){u});
}

#include "fractions.h"
#define cylinder (R - sqrt(sq(x - x0) + sq(y)))
event Cylinder (i++) {
  scalar cyl[];
  fraction (cyl, cylinder);
  foreach() {
    foreach_dimension()
      u.x[] -= cyl[]*u.x[];
    if (x < x0/2) {
      u.x[] = 1;
      u.y[] = 0;
    }
  }
  boundary ((scalar*){u});
}

event adapt (i++) 
  adapt_wavelet ((scalar*){u}, (double[]){0.02, 0.02}, maxlevel);

event mov (t += 0.5) {
  scalar omg[], cyl[];
  fraction (cyl, cylinder);
  foreach()
    cyl[] = 0.5 - cyl[];
  vorticity (u, omg);
  output_ppm (omg, file = "o.mp4", n = 600, box = {{x0/2, -6},{100, 6}},
	       min = -1.25, max = 1.25, linear = true, mask = cyl);
  foreach()
    omg[] = level;
  output_ppm (omg, file = "l.mp4", n = 600, min = 1, max = maxlevel,
	      box = {{x0/2, -5},{100, 5}});
}

event stop (t = 200);
