#include "grid/multigrid.h"
#include "diffusion.h"
#include "radial.h"

#include "run.h"
#include "view.h"
//#include "display.h"


scalar T[];
mgstats mgd;

const double Ri = 0.1; // inner radius
const double Ro = 1.;    // outer radius

int main() {
  size (Ro - Ri);
  origin (Ri, 0.);
  dtheta=2*pi;
  N=64;
  run();
}

T[right] = neumann(0.);
T[top] = neumann(0);

T[left] = dirichlet(0.5);
//T[bottom] = dirichlet(0.5); // useless ... don't have any impact

event init (i = 0) {
  foreach()
    T[] = 0.;
}

event diffusion (i++) {
  dt = dtnext(0.1);
  const face vector D[] = {0.01, 0.};
  mgd = diffusion (T, dt, D);
}

#define end 5
event profile (t = {0,end/2.,end}) {
  output_field ({T}, stdout, N);
}


void radial (coord * p) {
  double r = p->x, theta = p->y*dtheta/L0;
  p->x = r*cos(theta), p->y = r*sin(theta);
}

event video(i++){
  clear();
  view (quat = {0.000, 0.000, 0.000, 1.000}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = 0.004, ty = 0.023, tz = -4.074,
      width = 1280, height = 720);
  box();
  squares(color = "T", linear = "true");
  //cells();
  //isoline("T", n=21);
  save("T.mp4");
}

/**
## Results
![Temperature field.](hot-ring/T.mp4)
*/