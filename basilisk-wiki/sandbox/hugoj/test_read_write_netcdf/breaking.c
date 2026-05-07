/**

# Test reading/writing using netcdf_bas.h

This uses the breaking.c example
*/

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h" 
#include "layered/perfs.h"

//#include "view.h"
#include "bderembl/libs/netcdf_bas.h"

/**
The initial conditions are given by the wave steepness $ak$ and the
Reynolds number $Re=c\lambda/\nu$ with $c$ the phase speed of the
gravity wave and $\lambda$ its wavelength. */

double ak = 0.33;
double RE = 40000.;
double k_ = 2.*pi, h_ = 0.5, g_ = 1.;

#define T0  (k_*L0/sqrt(g_*k_))

int PAD=4; // padding for name of netcdf
int OUT=0; // index of the output file
char fileout[100];
/**
The domain is periodic in $x$ and resolved using 256$^2$
points and 30 layers. */

int main()
{
  origin (-L0/2., -L0/2.);
  periodic (right);
  N = 64;
  nl = 30;
  G = g_;
  nu = 1./RE;

  /**
  Some implicit damping is necessary to damp fast modes. This may be
  related to the slow/fast decoupling of the $\theta$-scheme mentioned
  by [Durran & Blossey, 2012](#durran2012). */
  
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
  /**
  We can use a larger CFL, in particular because we are not dealing
  with shallow-water/wetting/drying. */
  
  CFL = 0.8;

  /**
  The layer thicknesses follow a geometric progression, starting from
  a top layer with a thickness of 1/3 that of the regular
  distribution. */
  
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
  write_nc(); 
  OUT = OUT +1;

}

/**
We add (an approximation of) horizontal viscosity. */

event viscous_term (i++)
  horizontal_diffusion ((scalar *){u}, nu, dt);

/**
We log the evolution of the kinetic and potential energies.

~~~gnuplot Evolution of the kinetic, potential and total energy
set xlabel 't/T0'
plot [0:6]'log' u 1:2 w l t 'kinetic', '' u 1:3 w l t 'potential', \
     '' u 1:(($2+$3)/2.) w l t 'total/2'
~~~
*/

event logfile (i++; t<=T0) //t <= 8.*T0)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
    }
    gpe += sq(eta[])*dv();
  }
  //fprintf (stderr, "%g %g %g\n", t/T0, ke/2., g_*gpe/2.);
}

event end (i=end){
  write_nc(); 
  }




