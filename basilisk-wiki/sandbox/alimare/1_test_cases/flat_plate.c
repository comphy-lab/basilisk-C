/** 
![final_speed](flat_plate/ux.png)
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

#define MAXLEVEL 10
#define MINLEVEL 5

face vector muv[];

double u0 = 1.;
double h = 0.1;

double mygeom2(double x, double y){
  return y-h;
}

int main() {
  init_grid(256);
  size(1. [0]); // trick necessary for the use of embed ?
  mu = muv;
	run();
}

event init (t = 0) {
  CFL = 0.5;
  DT = 1.e-5;
  periodic(left);
	vertex scalar distn[];
	foreach_vertex() {
    distn[] = mygeom2(x,y);
  }
  boundary({distn});
  fractions(distn,cs,fs);
  fractions_cleanup(cs,fs);
  
  u.n[embed] = dirichlet(-u0); // u.n[embed] is actually x
}

event properties (i++) {
  foreach_face()
    muv.x[] = fs.x[]*u0;
}

event adapt (i++) {
	adapt_wavelet ({u,cs}, (double[]){1.e-2,1.e-2,1.e-2}, MAXLEVEL, MINLEVEL);
  fractions_cleanup(cs,fs);
}

event mylog(i+=10){
	stats s = statsf (u.x);
  stats s2 = statsf (u.y);
  fprintf (stderr, "%g %g %ld %g %g %g %g\n", \
  t, dt , grid->tn, s.min, s.max, s2.min, s2.max);
}

event final_event(i=100){
  view(tx = -0.5, ty = -0.5);
  draw_vof("cs","fs", filled = -1);
  squares("u.x");
  save("ux.png");
}
