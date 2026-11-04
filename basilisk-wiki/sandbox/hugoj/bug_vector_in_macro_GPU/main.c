double p0 = 1., rho = 1000.;

vector windmask[];

#define a_baro(f,k) (1.0*(-G*(f[k]-f[k-1])/Delta \
                - (p0/rho)*windmask.x[]*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))

#include "grid/cuda/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"


event init (t = 0)
{
  foreach() {
    zb[] = -2.;
    eta[] = 0.;
    double H = eta[] - zb[];
    foreach_layer()
      h[] = H*1/nl;
  }
}

event stop (t = 0.1) {}

int main()
{
  N = 16;
  nl = 4;
  L0 = 1.;
  origin (-L0/2.);
  periodic (right);
  G = 1.;
  run();
}
