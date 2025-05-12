/**
## test case radial diffusion
*/
#include "grid/multigrid.h"
#include "diffusion.h"
#include "radial.h"

#include "run.h"
#include "view.h"
//#include "display.h"


scalar T[];
mgstats mgd;

int main() {
  size (1.);
  dtheta=2*pi;
  N=64;
  run();
}

T[right] = neumann(0.);//adiabatic radial
T[top] = neumann(0);//??

event init (i = 0) {
  foreach()
    T[] = 1./sqrt(2*M_PI*sq(0.1))*exp(-sq(x)/(2*sq(0.1)));
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
  cells();
  isoline("T", n=21);
  save("T.mp4");
}


/**
## Results
~~~gnuplot
reset
set xlabel 'r'
set ylabel 'T'
p 'out' u ($2==0.0?$1:1/0):3 w p
~~~

![Temperature field.](radial-diff/T.mp4)

*/