/**
#Numerical setup of a half vortex ring through a stratified atmosphere

*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "navier-stokes/perfs.h"


scalar b[], * tracers = {b};

face vector av[];

double R = 1, Uo = 1;      // Normalized values;
double Re = 5000, Pi1 = 8, Pi4 = 0.01; // Dimensionless numbers
double nu, to, Nbv;             // Dependent variables
//double Nbv = 0.07;        // Statratification should depend on a new dimensionless number

#define RADIUS y

u.n[left] = dirichlet ((t <= to)*(RADIUS <= R)*Uo); // Injection
b[left] = dirichlet (0);

u.n[right] = t <= to ? neumann (0): dirichlet(0);
p[right]   = t <= to ? dirichlet (0): neumann (neumann_pressure(ghost));
b[right] = neumann (-sq(Nbv));
int maxlevel = 9;

int main() {
  L0  = 90*R;
  to  = Pi1*R/Uo;
  nu  = Uo*R/Re;
  Nbv = Pi4*Uo/R;
  a   = av;     // Link Gravity
  const face vector muc[] = {nu, nu};
  mu = muc;
  run();
}

event init (t = 0) {
  refine (RADIUS < R*2 && x < R/2 && level < maxlevel - 1); 
  refine (RADIUS < R   && x < R/5 && level < maxlevel);
  foreach()
    b[] = -sq(Nbv)*x;
  boundary ({b});
  DT = 0.1;
d}

event acceleration (i++) {
  boundary ({b});
  foreach_face(x) // x is the axial coordinate "z"
    av.x[] = -face_value(b,0);
  foreach_face(y) 
    av.y[] = 0;
}

event tracer_diffusion (i++) {
  diffusion (b, dt, mu);
}

event movies (t += 0.25) {
  scalar lev[], omg[];
  vorticity (u, omg);
  foreach()
    lev[] = level;
  output_ppm (omg, file = "omgpi40.01.mp4", n = 300, min = -1, max = 1);
  output_ppm (lev, file = "lpi40.01.mp4", n = 300, max = maxlevel);
  output_ppm (b, file = "bpi40.01.mp4", n = 300);
}

event track (t += 0.5) {
  scalar omg[];
  vorticity (u, omg);
  double xp = 0, yp = 0, tot = 0;
  foreach(reduction (+:xp) reduction (+:yp) reduction (+:tot)){
    if (omg[] > 0){
      xp += sq(omg[])*x*dv();
      yp += sq(omg[])*y*dv();
      tot += sq(omg[])*dv();
    }
  }
  if (tot != 0) {
    static FILE * fp = fopen ("datapi40.01", "w");
    fprintf (fp, "%g %g %g\n", t, xp/tot, yp/tot);
  }
}
  
event adapt (i++) {
  adapt_wavelet ({b, u}, (double[]){sq(Nbv)*R, 0.02, 0.02}, maxlevel );
}

event stop (t = 40*to);
