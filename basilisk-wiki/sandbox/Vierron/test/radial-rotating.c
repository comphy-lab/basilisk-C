#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "radial.h"
#include "navier-stokes/perfs.h"
#include "view.h"
//#include "display.h"

scalar f[];
scalar * tracers = {f};

face vector muc[];

const double Ri = 0.3; // inner radius
const double Ro = 1.;    // outer radius
double Re;

void radial (coord * p) {
  double r = p->x, theta = p->y*dtheta/L0;
  p->x = r*cos(theta), p->y = r*sin(theta);
}

int main() {
  size (Ro - Ri);
  origin (Ri, 0.);
  dtheta=2*pi;
  TOLERANCE=1e-3;
  DT = 1 [0];
  N=64;
  Re=1e2;
  mu=muc;
  run();
}


u.n[right] = dirichlet(0.);//u.r external boundary
u.t[right] = dirichlet(0.8);//u.theta external boundary

u.t[left] = dirichlet(1.);//u.theta internal boundary
u.n[left] = dirichlet(0.);//u.r internal boundary


event init (i = 0) {
  double angle1 = 20.*M_PI/180.;
  double r1 = (Ro-Ri)/4. + Ri;
  double r2 = 0.05;
  foreach()
    f[] = (sq(x*cos(theta) - (r1*cos(angle1))) + sq(x*sin(theta) -r1*sin(angle1)) < sq(r2)) ? 1. : 0;
}


event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*1./Re;
  boundary ((scalar*){muc});
}

#define end 10
event logfile (t = {0,end}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %f %g\n", t, s.sum, s.min, s.max);
}


#if 1
event video(i++){
  clear();
  view (quat = {0.000, 0.000, 0.000, 1.000}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = 0.004, ty = 0.023, tz = -4.074,
      width = 1280, height = 720);
  box();
  squares(color = "u.y", linear = "true");
  cells();
  save("uy.mp4");
}
#endif

/**
## Results
![uy field.](radial-rotating/uy.mp4)
*/