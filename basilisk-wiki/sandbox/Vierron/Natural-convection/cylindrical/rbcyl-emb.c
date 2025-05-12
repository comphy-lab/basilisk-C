/**
# Rayleight Benard instability in 3D with embed
*/
#include "grid/multigrid3D.h"
//#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "radial.h"
#include "navier-stokes/perfs.h"
#include "view.h"
//#include "display.h"


const int MAXLEVEL = 9;
const int MINLEVEL = 5;
const double Tbot = 0.5;
const double Ttop = -0.5;
const double Height = 1.;
const double EndTime = 40.;
int ii;

void radial (coord * p) {
  double r = p->x, theta = p->y*dtheta/L0;
  p->x = r*cos(theta), p->y = r*sin(theta);
}

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stderr, "%d %g %g %g %d\n", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
	    fflush(stderr);
}

scalar T[];
scalar * tracers = {T};
face vector muc[], av[];
mgstats mgT;

double Ra, Pr;

int main (int argc, char * argv[])
{
  size (55/15.);//R0
  dtheta=2*pi [0];
  Z0=-Height/2. + (L0/pow(2,MAXLEVEL+1))/2.;// to avoid superposition between grid and embedded boundary
  init_grid (1 << MINLEVEL);
  //unrefine(z>1.2*Z0 && level >= MINLEVEL);
  a=av;
  mu=muc;
  Ra = 1900 [0];
  Pr = 6.5 [0];
  run();
}

#if 1
//T[front] = dirichlet(-0.5);
T[back] = dirichlet(Tbot);
T[right] = neumann(0.);
T[embed] = dirichlet(Ttop);

u.n[embed] = dirichlet(0.);
u.n[right] = dirichlet(0.);
//u.n[front] = dirichlet(0.);
u.n[back] = dirichlet(0.);

u.t[embed] = dirichlet(0.);
u.t[right] = dirichlet(0.);
//u.t[front] = dirichlet(0.);
u.t[back] = dirichlet(0.);

u.r[embed] = dirichlet(0.);
u.r[right] = dirichlet(0.);
//u.r[front] = dirichlet(0.);
u.r[back] = dirichlet(0.);
#endif

event init (t=0) {
  solid (cs, fs, (+Height/2. + (L0/pow(2,MAXLEVEL+1))/2.)- z);
  foreach(){
    T[] = -z + (L0/pow(2,MAXLEVEL+1))/2.;
    foreach_dimension()
      u.x[] = 0.001*noise();
  }
  boundary ({T,u});
  DT = 0.001;
  dtnext(DT);
  TOLERANCE=10E-8;
   
  clear();
  view (quat = {-0.078, 0.599, 0.753, -0.261}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = -0.022, ty = 0.031, tz = -2.943,
      width = 1620, height = 960);
  box();
  squares(color = "T", n = {0,1,0}, alpha = 0.52, min=Ttop, max=Tbot);
  cells(n = {0,1,0}, alpha = 0.52);
  save("T.png");
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*Pr/sqrt(Ra);
  //boundary ((scalar*){muc});
}


event tracer_diffusion (i++) {
  //dt = dtnext (0.05*sq(1./pow(2,MINLEVEL))/(1./sqrt(Ra)));
  //dt = dtnext(0.005);
  face vector D[];
  foreach_face()
    D.x[] = fm.x[]*1./sqrt(Ra);
  //boundary ((scalar*){D});
  mgT = diffusion (T, dt, D);
  boundary ({T});
}

event acceleration (i++) {
  ii++;
  foreach_face(z)
    av.z[] += fm.z[]*Pr*(T[] + T[0,0,-1])/2.;
  foreach_face(y)
    av.y[] += fm.y[]*0.;
  foreach_face(x)
    av.x[] += fm.x[]*0.;
    
  if ((i==10||ii==10)){
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
  fprintf (stdout, "%f %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n", t, deltau, avg, du, s0.sum/s0.volume, s0.max, statsf(u.z).sum, normf(p).max, theta);
  fflush(stdout);
  mg_print (mgT);
  
}

#if 0
event adapt (i++){
  double err = 0.01;
  astats s = adapt_wavelet ((scalar*){cs, T, u},
		 (double[]){err*0.1, err, err, err, err}, MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
#endif

//numerical shadowgraph in cylindrical coordinate
event ombro (t = EndTime){
  scalar ombro[];
  double Deltax = 4./N;
  double Deltao = 2*M_PI/N;
  double deltaT = 1.;
  foreach(){
    ombro[] += (((T[1,0,0] - 2*T[] + T[-1,0,0])/sq(Deltax)) * deltaT + 1/x*(T[-1,0,0] - T[1,0,0])/(2*Deltax) + 1/sq(x)*(T[0,1,0] - 2*T[] + T[0,-1,0])/sq(Deltao) * deltaT) * (Deltax/sq(Height));
  }
view (quat = {0.000, 0.000, 0.000, 1.000}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = 0.094, ty = -0.045, tz = -3.934,
      width = 1920, height = 1080);
  squares(color = "ombro");
  save("ombro.png");
}

event movie (t += 1.; t<=EndTime) {
  clear();
  view (quat = {-0.078, 0.599, 0.753, -0.261}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = -0.022, ty = 0.031, tz = -2.943,
      width = 1920, height = 1080);
  squares (color = "pf", n = {0,0,1}, alpha = -0.45);
  box ();
  cells (n = {0,0,1}, alpha = -0.45);
  save("T.mp4");
}

event movie2(t += 1.; t<=EndTime){
  clear();
  view (quat = {-0.078, 0.599, 0.753, -0.261}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = -0.022, ty = 0.031, tz = -2.943,
      width = 1920, height = 1080);
  box ();
  cells ();
  isosurface (f = "T", color = "T", v=0);
  save("T2.mp4");
}

#if 0
#if TRACE > 1
event profiling (i += 19) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}

event profiling_plot (i = 40) {
  if (getenv ("DISPLAY"))
    popen ("gnuplot -e 'set term x11 noraise noenhanced title profiling' "
	   "$BASILISK/profiling.plot "
	   "& read dummy; kill $!", "w");
}
#endif
#endif

/**

## Results
![u.z field.](rbcyl-emb/T.mp4)
![isosurface T field.](rbcyl-emb/T2.mp4)

*/