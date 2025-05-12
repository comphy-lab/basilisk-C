/**
# Test the scale segment utility

![The length must be equal to the radius](scale/scale.png)
*/


#include "view.h"
#include "src/scale.h"

int main() {
  origin (-0.5, -0.5);
  init_grid (64);
  scalar f[];
  fraction (f, sq(x) + sq(y) - sq(0.2));
  view (fov = 30, tx = 0, ty = 0., width = 750, width = 750);
  draw_vof("f", filled = -1);
  box();
  scale(pos = {0., 0.3}, length = {0.2, 0.}, lw = 3, lc = {1, 0, 0});
  save("scale.png");
}
