/**
# A test for the particle advection with MPI

A test: there are about 100 times more particles than
cells.

![Many particles on 6 cores](test_manyparticles/mov.mp4)

*/
#include "navier-stokes/centered.h"
#include "view.h"
#define BVIEW 1
#include "particles.h"

int maxlevel = 10;
u.t[top] = dirichlet (0.);
#define RAD (sqrt(sq(x) + sq(y)))
#define ST (x/RAD)

const face vector nu[] = {0.005, 0.005};

int main() {
  mu = nu;
  L0 = 20.;
  X0 = Y0 = -L0/2.;
  N = 512;
  run();
}

event init (t = 0) {
  double k = 3.83170597;
  scalar psi[];
  foreach()
    psi[] = ((RAD > 1)*((1/RAD))*ST +
	     (RAD < 1)*((-2*j1(k*RAD)*ST/(k*j0(k))) + (RAD*ST)));
  boundary({psi});
  foreach() {
    u.x[] = -(psi[0, 1] - psi[0, -1])/(2*Delta);
    u.y[] = (psi[1] - psi[-1])/(2*Delta);
  }
  init_particles_in_cells();
}

event adapt (i++) {
  adapt_wavelet ((scalar*){u}, (double[]){0.05, 0.05}, maxlevel, 5);
  if (pid() == 0)
    printf ("%ld %d\n", grid->tn, sq(512));
}
  
event render (t += 0.5; t < 25) {
  scatter (loc, s = 5, pc = {sin(pid()), cos(pid()), sin(pid()*2.4)});
  box();
  save ("mov.mp4");
}
