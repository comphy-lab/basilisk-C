/**
# Tracers

First, `b` is the buoyancy and `c` is a tracer, next it is reversed. Finally, both a active tracers, but represent the buoyancy is different directions.

![`c`, `b` and particle tracers](ns4t/mov.mp4)
*/
#include "nsf4t.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

scalar b[], c[], * tracers = {b, c};
scalar zeros[]; // Just zeros

Particles parts;

int main() {
  foreach_dimension()
    periodic (left);
  L0 = 10;
  X0 = Y0 = -L0/2;
  const scalar muz[] = 0.001;
  kappa = nu = muz;
  N = 128;
  a.x = zeros;
  a.y = b;
  run();
  a.y = c;
  run();
  a.x = b;
  run();
}

event init (t = 0) {
  parts = init_tp_circle (100);
  foreach_particle_in(parts) {
    double R = sqrt(sq(x) + sq(y));
    p().x /= R;
    p().y /= R;
  }
  foreach_vert() {
    b[] = exp (-sq(x) - sq(y)) ;
    c[] = exp (-sq(x) - sq(y - 1))
      - exp (-sq(x) - sq(y + 1));
  }
}

event mov (t += 0.1) {
  view (width = 800, height = 300);
  squares ("b");
  translate (x = L0)
    scatter (parts);
  translate (x = -L0)
    squares ("c");
  save ("mov.mp4");
  
}


event stop (t = 15);

