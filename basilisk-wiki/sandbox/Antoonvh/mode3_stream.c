/**
# Mode 3 instability using the $\omega-\psi$ solver

![Mode 3 instability](mode3_stream/omg3.mp4)

See the [4th-order solver test](mode3f4.c)
 */

#include "navier-stokes/stream.h"

int maxlevel = 8;

int main() {
  foreach_dimension() {
    psi[left] = dirichlet(0);
    psi[right] = dirichlet(0);
  }
  L0 = 8.;
  origin (-pi*4./3., -2*exp(1) + 1.2); // Not centered
  N = 1 << maxlevel;
  run();
}

#define RAD (sqrt((sq(x) + sq(y))))
#define THETA(M) (M*asin(x/RAD))
#define RADP(P, M) ((1 + P*sin(THETA(M))))

double P = 1e-5, m = 3; // Perturbation and mode
double b = 1.45; //Vortex parameter

event init (t = 0) {
  refine (sq(x) + sq(y) < sq(b*1.2) && level < maxlevel + 2);
  double betav = 1./(1 - sq(b));//, alphav = -betav*sq(b);
  scalar omega_uf[];
  foreach() {
    double rp = RAD*RADP(P,m), omg = 0;
    if (rp <= 1.)
      omg = 2;
    else if (rp > 1 && rp <= b) 
      omg = 2*betav;
    omega_uf[] = omg;
  }
  //Helmholtz filter
  const face vector alphaf[] = {-sq(0.4/(2.*pi)), -sq(0.4/(2.*pi))};
  foreach_dimension()
    poisson (omega, omega_uf, alphaf, unity);  
}

event mov (t += 0.5) {
  output_ppm (omega, file = "omg3.mp4", map = cool_warm, n = 400);
}

event stop (t = 25);
