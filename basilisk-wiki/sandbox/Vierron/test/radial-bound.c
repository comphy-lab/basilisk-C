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
//T[top] = dirichlet(-0.5);// something weird happend

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

event boundarie (t=0){
  foreach_boundary(left){
    fprintf(stderr, "%s\n", "#left");
    fprintf(stderr, "%g %g %g %g %g\n", x, y, T[], T[-1,0], (T[] + T[-1,0])/2.);
    }
  fflush(stderr);  
  
  foreach_boundary(right){
    fprintf(stderr, "%s\n", "#right");
    fprintf(stderr, "%g %g %g %g %g\n", x, y, T[], T[1,0], (T[] + T[1,0])/2.);
    }
  fflush(stderr);  
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
  save("T.mp4");
}