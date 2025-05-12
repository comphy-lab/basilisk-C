#include "grid/multigrid3D.h"
#include "radial.h"
#include "run.h"
#include "view.h"
//#include "display.h"

int main() {
  size (npe()/2);
  dimensions(nx = npe()/2, ny = npe()/2, nz = npe()/4);
  dtheta=2*pi;
  N=64;
  run();
}

void radial (coord * p) {
  double r = p->x, theta = p->y*dtheta/L0;
  p->x = r*cos(theta), p->y = r*sin(theta);
}

event picture (t = 0) {
  clear();
  view (quat = {0.151, -0.258, -0.770, 0.564}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = -0.069, ty = 0.498, tz = -6.798,
      width = 1550, height = 940);
  cells ();
  box ();
  scalar procid[];
  foreach()
    procid[] = pid();
  squares("procid", min = 0, max = npe());
  save("cyl.png");
}
/**
![distribution pid()](cylmulti/cyl.png)
*/