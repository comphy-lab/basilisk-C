/**
# Rayleigh-Plateau instability

What happens when $R$ is varied? */

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

double R = 0.2, epsilon = 0.1;

int main() {
  f.sigma = 1.;
  rho2 = 0.01;
  N = 64;
  run();
}

event init (t = 0) {
  fraction (f, R*(1. + epsilon*cos(pi*x)) - y);
}

event gfsview (i += 10; t <= 10) {
  static FILE * fp = popen ("gfsview2D -s plateau1.gfv", "w");
  output_gfs (fp);
}
