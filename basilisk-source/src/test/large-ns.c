/**
# Large-amplitude standing wave

See [large.c]() for a more detailed description. */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"

int main()
{
  rho1 = 1./1000.;
  rho2 = 1.;
  G.y = -9.81;
  DT = 1e-3 [0,1];
  Y0 = - L0/2.;
  N = 256;
  run();
}

event init (i = 0)
{
  double a = 0.07, k = 2.*pi;
  fraction (f, - (a*cos(k*x) - y));
}

event profiles (t = 0.1; t += 0.1; t <= 0.5)
{
  output_facets (f, stderr);
}
