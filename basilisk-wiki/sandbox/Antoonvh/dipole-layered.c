/**
# Deformation of two vortex columns

Two vortex columns are initialized in a (rotating) fluid layer. Due
to surface drag, Ekman pumping (and suction) is generated. In turn,
the vortex collumns deform due to the in-plane divergence. 

Results without background rotation:

![Vorticity in the mid plane](dipole-layered/omg.mp4)

![Vertical velocity in the mid plane](dipole-layered/w.mp4)

Note: The vortex columns are unstable unless they are subject to
sufficient background rotation (see Taylor-Proudman theorem). 
*/
//#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"

double R = 0.15, U = 0.02, H = 0.3;
// size, speed and height
double Re = 500;
int main () {
  periodic (left);
  L0 =  10.*R;
  X0 = Y0 = -L0/2.;
  nl = 8;
  G = 10.;
  nu = U*R/Re;
  N = 128;
  run();
}

#define RAD (sqrt(sq(x) + sq(y)))
#define ST ((y)/RAD)

scalar psi[], omg[];
event init (t = 0) {
  double k = 3.8317/R;
  psi[top] = dirichlet (0.);
  psi[bottom] = dirichlet (0.);
  foreach() {
    omg[] = (RAD < R)*(-2*U*j1(k*RAD)*ST/(k*j0(k*R)))*sq(k);
    psi[] = 0;
  }
  poisson (psi, omg);
  boundary ({psi});
  foreach() {
    coord V, Vt = {0, 0};
    V.x = -(psi[0,1] - psi[0,-1])/(2*Delta);
    V.y = (psi[1] - psi[-1])/(2*Delta);
    vector u;
    scalar h;
    for ( h, u, in hl, ul) {
      h[] = H/nl;
      foreach_dimension()
	u.x[] = V.x + Vt.x;
    }
  }
#if TREE
  DT = 0.0025; //a small timestep for stability with adaptation
#endif
}

/**
The Coriolis force can be added by flipping a switch?? This creates an
asymetrical evolution.
*/
#if 0
double f = 1.;
event pressure (i++, last) {//Overload
  coord cor = {-f, f};
  vector u;
  face vector a;
  foreach_face() {
    for (a,u in al,ul) {
      a.x[] += cor.x*(u.y[] + u.y[-1])/2;
    }
  }
  boundary ((scalar*)al);
}
#endif
/**
The drag at the no-slip bottom is parameterized by an exponential decay law
with a timescale $\tau = \frac{h_0^2}{10*\nu}$.
 */
event bottom_drag (i++) {
  foreach() {
    vector u = ul[0];
    scalar h = hl[0];
    foreach_dimension()
      u.x[] *= exp(-dt*10*nu/sq(h[]));
  }
}

#if TREE
event adapt (i++) {
  vector u = ul[nl/2];
  boundary ({u.x, u.y});
  adapt_wavelet ({u.x, u.y}, (double[]){U/20, U/20}, 7);
}
#endif

event outputs (t += 0.1) {
  scalar omg[], omgb[];
  vorticity (ul[nl/2], omg);
  vorticity (ul[0], omgb);
  output_ppm (omg,  file = "omg.mp4",
	      min = -5*U/R, max = 5*U/R, n = 300);
  output_ppm (wl[nl/2], file = "w.mp4",
	      min = -U/10, max = U/10, n = 300);
  output_ppm (omgb, file = "omgb.mp4", min = -5*U/R, max = 5*U/R);
  output_ppm (eta,  file = "eta.mp4", min = 0.2999, max = 0.3001);
}

event stop (t = 20);
