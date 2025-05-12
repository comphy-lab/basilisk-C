#include "grid/multigrid.h"
#include "advection.h"
#include "radial.h"

#include "view.h"
//#include "display.h"

scalar T[];
scalar * tracers = {T};

void radial (coord * p) {
  double r = p->x, theta = p->y*dtheta/L0;
  p->x = r*cos(theta), p->y = r*sin(theta);
}



int main() {
  size (1.);
  dtheta=2*pi;
  // maximum timestep
  //DT = .1;
  N=64;
  run();
}


u.n[right] = dirichlet(0.);//u.r external boundary
u.t[right] = dirichlet(0.);//u.theta external boundary
//u.n[top] = dirichlet(0.);//useless
//u.t[top] = dirichlet(0.);//useless

T[right] = neumann(0.);

event init (i = 0) {
  foreach()
    T[] = (theta*180./M_PI > 0. && theta*180./M_PI < 20.) ? 1./sqrt(2*M_PI*sq(0.1))*exp(-sq(x-L0/2.)/(2*sq(0.1))) : 0;
}

event velocity (i++) {
  trash ({u});
  foreach_face(y)
    u.y[] = 0.01;
}

#define end 10.
event profile (t = {0,end/2.,end}) {
  output_field ({T}, stdout, N);
}


#if 1
event video(i++){
  clear();
  view (quat = {0.000, 0.000, 0.000, 1.000}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = 0.004, ty = 0.023, tz = -4.074,
      width = 1280, height = 720);
  box();
  squares(color = "T", linear = "true");
  cells();
  isoline("T", n=21);
  save("T.mp4");
}
#endif

/**
## Results
![Temperature field.](radial-advec/T.mp4)
*/