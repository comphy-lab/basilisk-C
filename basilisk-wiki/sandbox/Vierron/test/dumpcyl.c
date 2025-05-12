#include "grid/multigrid3D.h"
#include "radial.h"
#include "run.h"
#include "view.h"
//#include "display.h"

int main() {
  size (1.);
  Z0=-1./2;
  dtheta=2*pi;
  N=64;
  run();
}

event init(t=0){
  dump();
}

void radial (coord * p) {
  double r = p->x, theta = p->y*dtheta/L0;
  p->x = r*cos(theta), p->y = r*sin(theta);
}

event picture (t = 0) {
  view (quat = {0.486, 0.128, 0.237, 0.831}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = -0.014, ty = 0.088, tz = -2.988,
      width = 1550, height = 940);
  cells ();
  box ();
  squares (color = "cm");
  save("cyl.png");
}

/**
![Initial picture of cm](dumpcyl/cyl.png)
*/