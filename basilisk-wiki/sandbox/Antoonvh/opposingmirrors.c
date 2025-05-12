/**
![Gilles Deleuze demonstrates the opposing-mirror effect. Very heirarchical](https://i.pinimg.com/564x/c4/35/51/c435514f5ae4949e0d9f0ddcefa2ff63.jpg)

# Two opposing reflectors

![With a finite number of reflections](opposingmirrors/mirrors.png)

This is a test for the ray-depth counter in bwatch.
 */
#include "grid/multigrid3D.h"
#include "bwatch.h"

scalar s[];

int main() {
  X0 = Y0 = Z0 = -L0/2.;
  init_grid(N);
  foreach()
    s[] = x + y + z;
  boundary ({s});

  watch (fov = 1.5, O = {0.2, 0.3, 0.9}, poi = {-0.1, -0.1, -2});

  reflector (P = {0, 0, -1.});
  reflector (P = {0, 0, 1.},);
  dull_sphere(n = {0,0,1});
  sphere (s, C = {0, 0, -0.4}, R = 0.1, linear = true, map = gray);

  store (fopen ("mirrors.ppm", "w"));
  system ("convert mirrors.ppm mirrors.png");
}
