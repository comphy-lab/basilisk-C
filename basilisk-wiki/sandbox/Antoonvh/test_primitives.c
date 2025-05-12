/**
# Test the drawing of some primitives

![We draw an 3-color axis cross, a Gaussian, a cosine, an open-ended cylinder, and an arrow in the YZ plane](test_primitives/lines_etc.mp4)(loop)

*/
#include "grid/octree.h"
#include "view.h"
#include "some_primitives.h"

coord Gaussian (double ta) {
  return (coord){ta, exp(-sq(ta)), 5};
}

coord Cosinus (double ta) {
  return (coord){2.5, cos(2*ta), ta};
}

int main() {
  L0 = 10;
  init_grid (N);
  for (double theta = 0.1; theta < 2*pi + 0.1; theta += pi/60) { 
    view (theta = theta, phi = 0.3 + 0.2*sin(3*theta));
    box();
    draw_axis_cross();
    draw_paramterization (0.1, -5, 5, Gaussian, 100, (float[]){1., 0, 1.});
    draw_paramterization (0.1, -5, 5, Cosinus, 100, (float[]){.2, .4, 9.});
    draw_cylinder ({1,1,1}, {2,2,2}, end_cap = false);
    draw_arrow ({0, 1, 1}, {0, 3, 4}, len = 2, rad = 0.2, fc = {0.2, 0.2, 0.8});
    save ("lines_etc.mp4");
  }
}