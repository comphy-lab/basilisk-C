/**
# Cube of all colors

Our Eyes have [RGB
receptors](https://en.wikipedia.org/wiki/Cone_cell). We can visualize
the resulting color space with volumetric rendering.

![Do you see any new colors?](color_cube/cube.mp4)
 */
#include "grid/octree.h"
#include "bwatch.h"

int main() {
  X0 = Y0 = Z0 = -1.;
  L0 = 2;
  init_grid (32);
  scalar b[];
  scalar s[];
  vector colorv[];
  double ls = 1;
  double angh = 0.4;
  foreach() {
    s[] = exp((-sq(x) - sq(y) - sq(z))/ls);
    // Color code
    colorv.x[] = (x - X0)/L0;
    colorv.y[] = (y - Y0)/L0;
    colorv.z[] = (z - Z0)/L0;
  }
  multigrid_restriction({s, colorv});
  double fov = 5;
  double R = 10;
  FILE * fp = popen ("ppm2mp4 cube.mp4", "w");
  boundary ({s, colorv});
  for (double t = 0; t <= 2*pi; t += 2.*pi/(101)) {
    printf ("t %g\n", t);
    double theta = pi/2. + angh*sin(t);
    watch (O = {R*sin(theta)*cos(t), R*cos(theta), R*sin(theta)*sin(t)},
	   poi = {1e-6, 0, 0}, fov = fov);
    sphere (100, mat = {.dull = true});
    volume (s, cols = true, colorv = colorv, sc = .05, max = 1);
    store (fp);
    plain();
  }
  double t = 0;
  double theta = pi/2. + angh*sin(t);
  watch (O = {R*sin(theta)*cos(t), R*cos(theta), R*sin(theta)*sin(t)},
	 poi = {1e-6, 0, 0}, fov = fov);
  for (ls  = 1; ls >= 0.15; ls -= 0.025) {
    printf ("ls = %g\n", ls);
    foreach()
       s[] = exp((-sq(x) - sq(y) - sq(z))/ls);
    multigrid_restriction({s});
    sphere (100, mat = {.dull = true});
    volume (s, cols = true, colorv = colorv, sc = 0.05, max = 1);
    store (fp);
    plain();
  }
  for (double t = 0; t <= 2*pi; t += 2.*pi/(101)) {
    printf ("t %g\n", t);
    double theta = pi/2. + angh*sin(t);
    watch (O = {R*sin(theta)*cos(t), R*cos(theta), R*sin(theta)*sin(t)},
	   poi = {1e-6, 0, 0}, fov = fov);
    sphere (100, mat = {.dull = true});
    volume (s, cols = true, colorv = colorv, sc = 0.05, max = 1);
    store (fp);
    plain();
  }
}
