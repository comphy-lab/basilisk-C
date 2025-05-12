/**
# A dipole-wall collision with a simple NS solver

![Looks OK](tns/omg.mp4)

![Grid](tns/l.mp4)

~~~gnuplot
set xlabel 'Iteration'
set logscale y
set key top right
plot 'out' u 1:5 t 'Speed', '' u 1:6 t 'Cells'
~~~
*/
#include "ns.h"
#include "utils.h"

u.t[right] = dirichlet (0);
int maxlevel = 11;

int main () {
  periodic (bottom);
  L0 = 15;
  X0 = Y0 = -L0/2;
  nu = 0.0005;
  N = 1 << (maxlevel - 3);
  run();
}

#define RAD (sqrt(sq(x - xo) + sq(y - yo)))
#define ST ((y - yo)/RAD)

event init (t = 0) {
  CFL = 0.9;
  double k = 3.83170597, xo = 3, yo = 0;
#if TREE
  foreach_dimension()
    u.x.refine = refine_linear;
  refine (RAD < 2   && level < maxlevel - 1);
  refine (RAD < 1.1 && level < maxlevel);
#endif
  scalar psi[], omg[];
  foreach() 
    omg[] = -(RAD < 1)*-2*j1(k*RAD)*ST/(k*j0(k))*sq(k);
  poisson (psi, omg);
  foreach() {
    u.x[] = ((psi[0,1] - psi[0,-1])/(2*Delta));
    u.y[] = -(psi[1,0] - psi[-1,0])/(2*Delta);
  }
  boundary ({u.x, u.y});
}

event log_event (i += 10) 
  printf ("%d %g %d %d %g %ld\n", 
          i, t, mgp.i, mgp.nrelax, perf.speed, grid->tn);

event mov (t += 0.1) {
  scalar omg[], lev[];
  vorticity (u, omg);
  foreach() lev[] = level;
  output_ppm (omg, file = "omg.mp4",
	      min = -10, max = 10,     n = 512);
  output_ppm (lev, file = "l.mp4",
	      min = 1, max = maxlevel, n = 512);
}

#if TREE
event adapt (i++) 
  adapt_wavelet ((scalar*){u}, (double[]){0.01, 0.01}, maxlevel);
#endif

event stop (t = 15);
