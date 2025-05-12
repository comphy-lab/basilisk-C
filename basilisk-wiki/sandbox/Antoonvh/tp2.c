/**
# Test for `particles.h`

The result using 5 threads:

![Particles coloured with pid()-specific coding](tp2/locs.mp4)

The "shaking" is expected.

*/
vector u[];
#include "view.h"
#define BVIEW 1
#include "particles.h"

double xo, yo;
int main() {
  N = 8;
  run();
}

event init (t = 0) {
  P_RK3 = true;
  init_particles_in_cells();
  DT = 0.05;
}

event vel (i++) {
  xo = 0.5*sin(pi*t) + 0.55;
  yo = 0.5*cos(pi*t) + 0.35;
foreach() {
  u.x[] = (xo - x) ;
  u.y[] = (yo - y) ;
  }
 boundary ((scalar*){u});
}

#if BVIEW
event movie (i++) {
  view (tx = -0.5, ty = -0.5);
  vectors ("u", scale = 0.1);
  scatter (loc, pc = {sin(pid()), cos(pid()), sin(2.24*pid())});
  save ("locs.mp4");
}
#endif

event set_dt (i++) {
  dt = dtnext(DT);
}
  
event stop (t = 10);


