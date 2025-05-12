/**
# Vortex ring injection into a stratified atmosphere 

![Buoyancy](axiringstrat/b.mp4)

Warm air is transported along the axial coordinate, but the ring
dissapears as a buoyant tail exhaust warm fluid.

![vorticity](axiringstrat/omg.mp4)

## todo

* Introduce a dimensionless number that compares the stratification
 ($N\ [1/s]$) against vortex-injection parameters $U$ and $R$  
* Add code to identify vortex position  
* Analysis  
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar b[], * tracers = {b};

face vector av[];

double R = 1, Uo = 1;      // Normalized values;
double Re = 5000, Pi1 = 8; // Dimensionless numbers
double nu, to;             // Dependent variables
double Nbv = 0.07;        // Statratification should depend on a new dimensionless number

#define RADIUS y

u.n[left] = dirichlet ((t <= to)*(RADIUS <= R)*Uo); // Injection
b[left] = dirichlet (0);

u.n[right] = t <= to ? neumann (0): dirichlet(0);
p[right]   = t <= to ? dirichlet (0): neumann (neumann_pressure(ghost));
b[right] = neumann (-sq(Nbv));
int maxlevel = 9;

int main() {
  L0  = 25*R;
  to  = Pi1*R/Uo;
  nu  = Uo*R/Re; 
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
}

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
  output_ppm (omg, file = "omg.mp4", n = 300, min = -1, max = 1);
  output_ppm (lev, file = "l.mp4", n = 300, max = maxlevel);
  output_ppm (b, file = "b.mp4", n = 300);
}

event adapt (i++) {
  adapt_wavelet ({b, u}, (double[]){sq(Nbv)*R, 0.02, 0.02}, maxlevel );
}

event stop (t = 6*to);
