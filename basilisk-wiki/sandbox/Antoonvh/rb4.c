/**
# 2D Rayleigh-Bernard convection in a periodic square. 

Setting $Ra = 10^6$ and $Pr = 1$ results in a periodic-in-time
pattern.

![Buoyancy](rb4/b.mp4)
 */
#define NOSLIP_TOP (1)
#define NOSLIP_BOTTOM (1)
#include "nsf4t.h"
scalar zeros[], b[], * tracers = {b};

b[bottom] = dirichlet_vert_bottom(0.5);
b[top]    = dirichlet_vert_top(-0.5);

int main() {
  periodic (left);
  const scalar diff[] = 1e-3;
  kappa = nu = diff;
  a.x = zeros;
  a.y = b;
  run();
}
 
event init (t = 0) {
  TOLERANCE = 1e-5;
  foreach_vert()
    b[] = noise()/100;
  boundary ({b});
}

event mov (t += 0.1) 
  output_ppm (b, file = "b.mp4", n = 300, min = -0.55, max = 0.55);

event stop (t = 100);
