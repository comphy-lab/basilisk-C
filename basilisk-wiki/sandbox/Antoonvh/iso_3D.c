/**
# Test isolines in 3D

![This is the intended result](iso_3D/lines.mp4)

Cross section of a volume fraction field
![This is the intended result](iso_3D/cross.mp4)
*/
#include "grid/octree.h"
#include "view.h"
#include "iso3D.h"

scalar s[];

int main() {
  X0 = Y0 = Z0 = -L0/2;
  init_grid (16);
  refine (x < 0 && level < 5);
  unrefine (z < 0 && y < 0 && level > 2);
  /**
  The function is tested under a bunch of angles
  */
  for (double theta = 0; theta <= 3*pi; theta += 0.02) {
    foreach()
      s[] = exp(-sq(x* (1 + 0.05*sin(theta)))
		-sq(y* (1 + 0.2 *cos(theta)))
		-sq(z/ (1 + 0.1 *cos(theta))));
    boundary ({s});
    coord np = {sin(theta), cos(0.5*theta), sin(1.5*theta)};
    double alphap = 0.155435 + 0.2*sin(2*theta);
    view (theta = theta/10, phi = sin(theta)/4., psi = cos(theta)/8.);
    isoline2 ("s", val = 0.75 + 0.05*sin(theta),
	      np = np, alpha = alphap, lc = {1,0,1}, lw = 3);
    isosurface ("s", 0.78 + 0.05*sin(theta));
    box();
    cells (n = np, alpha = alphap);
    save ("lines.mp4");
  }
  scalar f[];
  fraction (f, sq(1.2*x) + sq(y/1.2) + sq(z) - sq(0.3));
  for (double theta = 0 ; theta < 2*pi; theta += 0.02) { 
    view (theta = theta, phi = 0.01, psi = 0.01);
    cross_section ("f", alpha = 0.05);
    save ("cross.mp4");
  }
}
