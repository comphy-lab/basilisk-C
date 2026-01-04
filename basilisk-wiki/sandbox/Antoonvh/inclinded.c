/**
# Intermittent chaos in a stratified inclinded

![](inclinded/mov.mp4)

![The grid looks fine](inclinded/cells.mp4)
 */
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"

// Only works in a centered domain
void tube_geometry (double L, double H, double rad) {
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = H/2. + rad - fabs(y); // duct ...
    phi[] = max (phi[], fabs(x) - L/2. + rad); // ... or tank
    if (fabs(y) > H/2 && fabs(y) < (H/2 + rad) && fabs(x) > L/2 - rad && fabs(x) < L/2)
      phi[] = sqrt(sq(fabs(y) - (H/2 + rad)) + sq(fabs(x) - (L/2 - rad))) ;
  }
  fractions (phi, cs, fs, rad);
  output_ppm (phi, file = "phi.png", min = -rad, max = 2*rad);
}

//buoyancy field
scalar b[], * tracers = {b};
//angle of inclination
double incline = pi*5./180.;

face vector av[], muc[];

u.n[embed] = dirichlet (0);
u.t[embed] = dirichlet (0);

double nu = 1/2000.; // = kappa

int maxlevel = 9;
// Geometry
double Le = 6, H = 0.5, rad = 0.25;

int main() {
  L0 = 10;
  X0 = Y0 = -L0/2;
  N = 1 << 6;
  a = av;
  mu = muc;
  DT = 0.1;
  run();
}
/**
Fluid properties are set in the fluid domain
 */
event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*nu;
}

double k = 5, b0 = 1;

#define STRAT(x) (-b0*tanh(k*x)*(cs[] > 0))

event init (t = 0) {
  /**
     Initialize the domain and stratification.
   */
  astats as = {.nf = 9999, .nc = 9999};
  while (as.nf > 20 || as.nc > 20) {
    tube_geometry(Le, H, rad);
    foreach()
      b[] = (-b0*tanh(k*x)*(cs[] > 0));
    as = adapt_wavelet({cs, b, u}, (double[]){1e-2, 5e-2, 2e-2, 2e-2}, maxlevel, 5);
  }
  // Slope limiter
  b.gradient = minmod2;
}

event acceleration (i++) { //Buoyancy
  coord dir = {sin(incline), cos(incline)};
  foreach_face()
    av.x[] = dir.x*(b[] + b[-1])/2.;
}

event tracer_diffusion (i++) {
  diffusion (b, dt, mu);
}
/**
Dampen the flow in the tanks
*/
event dampen (i++, last) {
  double tau = .5;
  foreach() {
    if (fabs(x) > Le/2 + H || fabs(y) > (H)) {
      b[] -= dt/tau*(b[] - -b0*tanh(k*x)*(cs[] > 0));
      foreach_dimension()
	u.x[] -= 0.5*dt/tau*u.x[];
    }
  }
}
event movie (t += 0.5) {
  output_ppm (b, file = "b.mp4", n = 300);
  view (psi = -incline);
  draw_vof ("cs", "fs", filled = -1, fc = {0.6, 0.6, 0.6});
  draw_vof ("cs", "fs", lw = 2);
  squares ("b", min = -1.1, max = 1.1);
  save ("mov.mp4");
  save ("snapshot.png");
  cells();
  save ("cells.mp4");
  save ("cells.png");

}

event adapt (i++) {
  adapt_wavelet({cs, b, u}, (double[]){1e-2, 5e-2, 2e-2, 2e-2}, maxlevel, 5);
}

event stop (t = 100);

