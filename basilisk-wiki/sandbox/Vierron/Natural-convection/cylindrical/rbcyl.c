/**
# Rayleight Benard instability in 3D

don't work with adaptative mesh
*/

//#include "grid/octree.h"
#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "radial.h"
#include "navier-stokes/perfs.h"
#include "view.h"
//#include "display.h"

#define MAXLEVEL 8
#define MINLEVEL 5

scalar T[];
scalar * tracers = {T};

face vector muc[], av[];

double Ra, Pr;

int main (int argc, char * argv[])
{
  size (4.);
  dtheta=2*pi;
  Z0=-L0/2.;
  init_grid (1 << 4);
  a=av;
  mu=muc;
  Ra = 1900; Pr = 0.71;
  run();
}

T[front] = dirichlet(-0.5);
T[back] = dirichlet(0.5);
T[left] = neumann(0.);
T[right] = neumann(0.);
T[top] = neumann(0.);
T[bottom] = neumann(0.);


u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.n[front] = dirichlet(0.);
u.n[back] = dirichlet(0.);

u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.t[front] = dirichlet(0.);
u.t[back] = dirichlet(0.);

u.r[top] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);
u.r[right] = dirichlet(0.);
u.r[left] = dirichlet(0.);
u.r[front] = dirichlet(0.);
u.r[back] = dirichlet(0.);

event init (t=0) {
  foreach(){
    T[] = -z/(L0/2.);
    foreach_dimension()
      u.x[] = 0.01*noise();
  }
  boundary ({T,u});
  //dump("init");
  DT = 0.001;
  dtnext(DT);
  TOLERANCE=10E-6;
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*Pr/sqrt(Ra);
  boundary ((scalar*){muc});
}


event tracer_diffusion (i++) {
  face vector D[];
  foreach_face()
    D.x[] = fm.x[]*1./sqrt(Ra);
  boundary ((scalar*){D});
  diffusion (T, dt, D);
  boundary ({T});
}

event acceleration (i++) {
  foreach_face(z)
    av.z[] += fm.z[]*Pr*(T[] + T[0,0,-1])/2.;
  foreach_face(y)
    av.y[] += 0.;
  foreach_face(x)
    av.x[] += 0.;
  if ((i==10||i==10)){
    DT=0.1;
    TOLERANCE=10E-4;
  }
}

scalar un[];
event init_un (i = 0) {
  foreach()
    un[] = u.z[];
}

event logfile(i++){
  scalar div[];
  double deltau = change (u.z, un);
  double avg = normf(u.z).avg, du = deltau/(avg + SEPS);
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += u.x[1] - u.x[];
    div[] /= Delta;
   }
  stats s0 = statsf (div);
  fprintf (stdout, "%f %.9g %.9g %.9g %.9g %.9g %.9g %.9g \n", t, deltau, avg, du, s0.sum/s0.volume, s0.max, statsf(u.z).sum, normf(p).max);
  fflush(stdout);
}

#if 0
event adapt (i++){
  double err = 0.01;
  astats s = adapt_wavelet ((scalar*){T, u},
		 (double[]){err, err, err, err}, MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
#endif

void radial (coord * p) {
  double r = p->x, theta = p->y*dtheta/L0;
  p->x = r*cos(theta), p->y = r*sin(theta);
}

event movie (t += 1.; t<=40.) {
  clear();
    view (quat = {-0.206, 0.558, 0.759, -0.266}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = 0.115, ty = -0.025, tz = -3.580,
      width = 1548, height = 936);
  box ();
  cells ();
  //isosurface (f = "T", color = "T", min = -0.5, max = 0.5);
  squares (color = "T", min = -0.5, max = 0.5, n = {1,0,0}, alpha = 2);
  cells (n = {1,0,0}, alpha = 2);
  save ("T.mp4");
}

event movie2(t += 1.; t<=40.){
  clear();
  view (quat = {0.183, 0.612, 0.731, 0.242}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = 0.109, ty = 0.152, tz = -2.879,
      width = 1980, height = 1020);
  box ();
  isosurface (f = "T", v = -0.45, min = -0.5, max = 0.5);
  isosurface (f = "T", v = 0.45, min = -0.5, max = 0.5);
  isosurface (f = "T", v = 0., min = -0.5, max = 0.5);
  save ("T2.mp4");
}

event movie3(t+=1.; t<=40.){
  clear();
  view (quat = {0.544, 0.239, 0.173, 0.785}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = -0.104, ty = -0.054, tz = -2.922,
      width = 1548, height = 936);
  squares (color = "u.z", max = 0.9, spread = 1);
  box ();
  cells ();
  squares (color = "T", min = -0.48, max = 0.5, n = {1,0,0}, alpha = 1);
  save ("T3.mp4");
}

/**
## Results
![video1](rbcyl/T.mp4)
![video2](rbcyl/T2.mp4)
![video3](rbcyl/T3.mp4)
*/