/**
# 3D breaking wave

*/

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"

#if OUTNC
#include "bderembl/libs/netcdf_bas.h"
#endif


double ak = 0.5;
double RE = 20000.;
double k_ = 2.*pi, h_ = 1., g_ = 1.;

#include "hugoj/lib/common_waves.h"

#define T0  (2.*pi/sqrt(g_*k_))
double lam = 0. [1];          // Wave length
double NT0 = 4. [0, 0, 0];   // Duration of the experiment
double cp = 0.[1, -2];        
double omega = 0. [0,-1];

/**
The domain is periodic in $x$ and resolved using 256$^2$
points and 30 layers. */

int main()
{
  origin (-L0/2., -L0/2.);
  periodic (right);
  N = 64;
  nl = 64;
  G = g_;

  // Wave
  lam = 2*pi/k_;
  omega = sqrt(g_*k_);
  cp = omega/k_;
  nu = cp*lam/RE;
  DT=0.005;

  /**
  Some implicit damping is necessary to damp fast modes. This may be
  related to the slow/fast decoupling of the $\theta$-scheme mentioned
  by [Durran & Blossey, 2012](#durran2012). */
  
  fprintf(stderr, 
      "T0 = %f, lam=%f, g_=%f,omega=%f,nu=%g\n", T0, lam, g_, omega, nu);
  run();
}


event init (i = 0)
{
  geometric_beta (1./3., true);

  foreach() {
    zb[] = - h_;
    double H = wave_stokes(x, 0, ak, k_, h_) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      z += h[]/2.;
      u.x[] = u_x_stokes(x, z, ak, k_, h_);
      w[] = u_y_stokes(x, z, ak, k_, h_);
      z += h[]/2.;
    }
  }
  #if OUTNC
  create_nc({zb, eta, h, u, w}, "out.nc");
  #endif
}

/**
We add (an approximation of) horizontal viscosity. */

event viscous_term (i++) {
  horizontal_diffusion ({u, w}, nu, dt);
}

/**
We log the evolution of the kinetic and potential energies.

~~~gnuplot Evolution of the kinetic, potential and total energy
set xlabel 't/T0'
plot [0:6]'log' u 1:2 w l t 'kinetic', '' u 1:3 w l t 'potential', \
     '' u 1:(($2+$3)/2.) w l t 'total/2'
~~~
*/

event logfile (i++; t <= NT0*T0)
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
  fprintf (stderr, "%g %g %g\n", t/T0, ke/2., g_*gpe/2.);
}

/**
And generate the movie of the free surface (this is quite
expensive). The movie is 45 seconds at 25 frames/second. 
 */

#if MOVIE
event movie (t += NT0*T0/(45*25))
{
  vector u = lookup_vector ("u29");
  output_ppm (u.x, file = "movie.mp4", linear = true, min = -0.15, max = 0.6, n = 512);
}
#endif 

/**
  Optional: output fields for analysis
*/
#if OUTNC
event writenc (i++; t <= NT0*T0) 
{
  if (i%50==0){
    write_nc();
  }
}
#endif

