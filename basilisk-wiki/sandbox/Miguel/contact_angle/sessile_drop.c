/**
Liquid droplet on a solid substrate. The wettability of the substrate can be changed via theta_0
*/
#define theta_0 (pi/4.)
#include "navier-stokes/centered.h"
#include "two-phase_angle.h"
#include "tension_angle.h"
#include "axi.h"

 double R = 1;
 vector h[];
 u.t[left] = dirichlet(0.);

int main()
{
  init_grid (32);
  origin (-0.0, -0.0);
  mu1 = 1.; 
  mu2 = 1.; 
  rho1 = 2.; 
  rho2 = 1.; 
  L0 = 4*R; 
  run();

}

event init (t=0){
  f.sigma = 32.;
  fraction (f, R - sq(x) - sq(y));
}

event profile(t = 0; t < 2; t+= 5e-2) {
  char name[80];
  sprintf(name,"profile_t%f.dat", t);
  FILE * fp1 = fopen(name,"w");
  output_facets (f, fp1);
}

/**
## Results

Droplet shape as the time advances:

~~~gnuplot Facets of the droplets
set size ratio -1
plot 'profile_t0.000000.dat' w l, \
'profile_t0.500000.dat' w l, \
'profile_t1.000000.dat' w l, \
'profile_t1.500000.dat' w l, \
~~~
 */
