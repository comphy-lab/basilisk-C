/**

# Test reading/writing using netcdf_bas.h

This uses the breaking.c example
*/

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h" 
#include "layered/perfs.h"

#include "bderembl/libs/netcdf_bas.h"

double ak = 0.33;
double RE = 40000.;
double k_ = 2.*pi, h_ = 0.5, g_ = 1.;

double T0=1.0; //(k_*L0/sqrt(g_*k_))

int PAD=4; // padding for name of netcdf
int OUT=0; // index of the output file
char fileout[100];


int main()
{
  origin (-L0/2., -L0/2.);
  periodic (right);
  N = 64;
  nl = 30;
  G = g_;
  nu = 1./RE;
  theta_H = 0.51;

  run();
}

/**
The initial conditions for the free-surface and velocity are given by
the third-order Stokes solution. */

#include "test/stokes.h"

event init (i = 0)
{
  fprintf(stderr, "T0=%f\n",T0);
  CFL = 0.8;

  geometric_beta (1./3., true);

  double a = 0.25;
  foreach() {
    zb[] = - 0.5 + a*sin(k_/2.*y);
    double H = wave(x, 0) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      z += h[]/2.;
      u.x[] = u_x(x, z);
      w[] = u_y(x, z);
      z += h[]/2.;
    }
  }
  
  // Writing ini step
  sprintf(fileout, "%0*d.nc", PAD, OUT); // add the padding
  create_nc({zb, h, u, w}, fileout);
  fprintf(stderr, "after create_nc\n");
  write_nc(); 
  fprintf(stderr, "after write_nc\n");
  OUT = OUT +1;
}

event end (t=T0){
  write_nc(); 
  }




