/**
[<img src="https://climateextremes.org.au/wp-content/uploads/2020/05/Underwater-by-Alexandra-Kock-Pixabay-1280x640.jpg" alt="drawing" width=80%/>](https://climateextremes.org.au/research-brief-why-ocean-models-fail-to-replicate-turbulent-convection/)  
Vortices can entrain air bubbles (click image for link to the climate extremes website)

# Vortex-air entrainment

![Entrainment occurs](lamb-air/mov.mp4)
 */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "view.h"

#define RAD (sqrt(sq(x) + sq(y)))
#define ST (x/RAD)

double h = 5, grav = 2;

int maxlevel = 9;
double ue = 1e-3;

int main() {
  L0 = 20;
  X0 = Y0 = -L0/2.;
  const face vector G[] = {0, -grav};
  mu2 = 1e-3;
  rho2 = 1;
  mu1 = mu2;
  rho1 = 10;
  a = G;
  run();
}

event init (t = 0) {
  TOLERANCE = 1e-4;
  fraction (f, h - y);
  refine (RAD < 2 && level < maxlevel);
  double k = 3.83170597;
  scalar psi[];
  foreach()
    psi[] = ((RAD > 1)*((1/RAD))*ST +
	     (RAD < 1)*((-2*j1(k*RAD)*ST/(k*j0(k))) + (RAD*ST)));
  boundary({psi});
  coord f = {-1.,1.};
  foreach() {
    foreach_dimension()
      u.x[] = f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
  }
  boundary ((scalar*){u});
}

event mov (t += 0.1; t < 25) {
  scalar omg[];
  vorticity (u, omg);
  foreach() 
    omg[] = fabs(omg[]) < 5e-1 ? nodata : omg[];
  squares ("omg", min = -5, max = 5);
  draw_vof ("f", fc = {0.3, 0.3, 0.8}, filled = 1);
  draw_vof ("f", lw = 2);
  save ("mov.mp4");
}

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.001, ue, ue}, 8, 6);
}
