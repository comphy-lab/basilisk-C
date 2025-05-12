/**
# The big bang

Mass fragements are given an initial spin, leading to the formation
of a few galaxies.

![A movie of the mass elements](bigbang/galaxies.mp4)
 */
#include "grid/octree.h"
#include "cosmology.h"
#include "view.h"
#include "scatter2.h"

double Omg = 2;

int main() {
  L0 = 100;
  X0 = Y0 = Z0 = -L0/2;
  DT = 0.01;
  run();
}

event init (t = 0) { 
  cosmos = init_ip_square (40, l = 7);
  foreach_particle_in(cosmos) {
    p().m = 1 + 0.1*noise();
    p().u.x =  Omg*y;
    p().u.y = -Omg*x;
    p().u.z = 0;
    p().z   = 1 + noise();
  }
}

event set_dtmax (i++)
  dt = dtnext(DT);

event movie (t += 0.05) {
  view (fov = 8, theta = 0.5*sin(t/4), phi = 0.2*cos(t/3));
  scatter (cosmos);
  save ("galaxies.mp4");
}

event adapt (i++)
  adapt_wavelet ({Phi}, (double[]){0.2}, 8);

event stop (t = 10);
