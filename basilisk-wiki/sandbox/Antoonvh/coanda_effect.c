/**
# Coandă effect

![A jet may bend along an object. The movie shows an object (gray), the velocity-field magnitude (white to red), and tracer particles](coanda_effect/mov.mp4)
*/
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "tracer-particles.h"
#include "scatter2.h"

Particles tracer_p;
// Jet
u.n[left] = dirichlet (fabs(y) < 0.5);
u.t[left] = dirichlet (0.);
p[left] = neumann (0);
// Outflow
u.n[right] = neumann (0);
p[right] = dirichlet (0);
// Object
u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);

face vector muc[];
double muv = 1./100.;
int main() {
  L0 = 15;
  X0 = Y0 = -L0/2.;
  N = 256;
  mu = muc;
  run();
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*muv;
}

event init (t = 0) {
  // particles
  tracer_p = new_tracer_particles(0);
  // Geometry
  vertex scalar phi[];
  coord corner = {-6, -5.75};
  double radius = 5;
  foreach_vertex() {
    if (x < corner.x)
      phi[] = y - corner.y;
    else if (y < corner.y)
      phi[] = x - corner.x;
    else
      phi[] = sqrt(sq(x - corner.x) + sq(y - corner.y));
  }
  fractions (phi, cs, fs, radius);
  // Initial jet
  foreach()
    if (fabs(y) < 0.5)
      u.x[] = 1;
}
/**
Particles are periodically added and removed once they approach the outlet boundary
 */
event add_particles (t += 2) {
  particle p1 = {0};
  p1.x = X0;
  p1.y = noise()/2;
  set_a_particle_attributes(&p1);
  add_particle (p1, tracer_p);
}

event and_remove_particles (i++) {
  remove_particles (tracer_p, x > X0 + L0*0.99);
}

event movies (t += 0.5) {
  scalar U[];
  draw_vof ("cs", "fs", fc = {0.3, 0.3, 0.3}, filled = -1);
  
  foreach()
    U[] = sqrt(sq(u.x[]) + sq(u.y[]));
  squares ("U", map = blue_white_red, linear = true, min = -0.75, max = 0.75);
  draw_vof ("cs", "fs", lw = 2);
  scatter (tracer_p);
  save ("mov.mp4");
  
}
double ue = 0.02;
event adapt (i++) {
  adapt_wavelet({cs, u}, (double[]){1e-3, ue, ue}, 8);
}

event stop (t = 250);
