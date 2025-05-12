/**
# A soil underneath an atmosphere

On this page we explore a possibility to implement a soil underneath
an atmosphere. The soil may store and diffuse heat with a different
diffusivity compared to the air in the atmosphere.

The case concerns warm air that is transpored downwards by a dipolar
vortex from aloft. The warm air reaches the underlaying surface and the
heat then partly diffuses into the soil where it is stored. 

![This is the movie we generate](soil/mp4.mp4)

It *appears* to work OK.
*/
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "fractions.h"
#include "view.h"

#define RAD (pow(pow((x - xo), 2)+(pow((y - yo), 2)), 0.5))
#define ST (-(x - xo)/RAD)

scalar s[], f[];
scalar * tracers = {s};
s[top] = neumann(1.);
double xo = 7.7, yo = 7.6;
double temp = 30;
double Re = 1000;

int main(){
  L0 = 15.;
  init_grid (1 << (8));
  f.refine = f.prolongation = fraction_refine;
  foreach_dimension()
    u.x.refine = refine_linear;
  const face vector muc[] = {1./Re, 1./Re};
  mu = muc;
  run();
}

event init (t = 0){
  fraction(f,  pi - y);
  foreach()
    s[] = (y - pi)*(y > pi);
  refine (RAD < 2.0 && level <= 8);
  refine (RAD < 1.0 && level <= 9);
  scalar psi[];
  double k = 3.83170597;
  foreach() 
    psi[] = ((RAD > 1)*((1/RAD))*ST) + ((RAD < 1)*((-2*j1(k*RAD)*ST/(k*j0(k))) + (RAD*ST)));
  boundary({psi});
  foreach() {
    u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta));
    u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta);
  }
  boundary(all);
}

face vector muz[];
event tracer_diffusion(i++){
  foreach_face()
    muz.x[] = (1./Re) + 0.05*(f[] + f[-1]);
  diffusion(s, dt , muz);
}

event adapt(i++)
  adapt_wavelet((scalar *){u, s}, (double []){0.05, 0.05, 0.025}, 9); 

event wall(i++){ //Stephane's trick
  foreach(){
    foreach_dimension()
      u.x[] -= u.x[]*f[];
  }
}

event bviewer(t += 0.075; t <= temp){
  scalar omega[];
  view(fov = 25, tx = -0.5, ty = -0.4, width = 1200, height = 500);
  vorticity(u, omega);
  squares("omega", map = cool_warm);
  draw_vof("f", lw = 3);
  translate(x = L0){
    squares("s", min = 0, max = 1);
    draw_vof("f", lw = 3);
  }
  translate(x = -L0){
    cells();
  }
  draw_string(" Cells, Vorticity and Temperature", size = 35);
  save("mp4.mp4");
}
