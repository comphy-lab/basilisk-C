/**
# A dipole-wall collision with a simple staggered NS solver on an adaptive grid

![Looks OK](tnsfa/omg.mp4)

![Grid level](tnsfa/lev.mp4)
 */
#include "nsf.h"

u.t[top] = dirichlet (0);

int main() {
  periodic (left);
  L0 = 15;
  X0 = Y0 = -L0/2;
  nu = 0.0002;
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
/**
For adaptation we compute a centered velocity field.
 */
event adapt (i++) {
  vector uc[];
  uc.t[top] = dirichlet (0.);
  boundary ((scalar*){u});
  foreach()
    foreach_dimension()
      uc.x[] = (u.x[] + u.x[1])/2.;
  boundary ((scalar*){uc});
  adapt_wavelet ({uc.x, uc.y}, (double[]){1e-2, 1e-2}, 10);
}

event mov (t += 1) {
  scalar omg[], lev[];
  vorticity (u, omg); //never mind...
  foreach()
    lev[] = level;
  output_ppm (omg, file = "omg.mp4");
  output_ppm (lev, file = "lev.mp4");
}

event stop (t = 200);
