/**
# A dipole-wall collision with a simple staggered NS solver

![Looks OK](tnsf/omg.mp4)

 */
#include "grid/multigrid.h"
#include "nsf.h"

u.t[top] = dirichlet (0);

int main () {
  periodic (left);
  L0 = 15;
  X0 = Y0 = -L0/2;
  nu = 0.002;
  N = 256;
  run();
}

void add_vortex (double x0, double y0, double s) {
 foreach_face(x) 
   u.x[] +=  s*(y - y0)*exp(-sqrt(sq(x - x0) + sq(y - y0)));
 foreach_face(y)
   u.y[] += -s*(x - x0)*exp(-sqrt(sq(x - x0) + sq(y - y0)));
}

event init (t = 0) {
  add_vortex ( 1, 0,  1);
  add_vortex (-1, 0, -1);
}

event mov (t += 1) {
  scalar omg[];
  vorticity (u, omg);
  output_ppm (omg, file = "omg.mp4");
}

event stop (t = 200);
