/**
# Collision of 2 droplets in a periodic box

This is a case to test the three-phase Navier-Stokes solver. Results should be
the same as the case with two-phases and [no-coalescence.h](no-coalescence.h).

This case is made non-dimensional using the drop radius $R$ and drop velocity
$U$. The problem parameters are
$$
\mathrm{Re} = 1460, \quad 
\mathrm{We} = 0.1, \quad 
\rho^* = 833, \quad 
\mu^* = 55 
$$

Here, colors indicate the volume fractions.

![Results from `collision2.c`](collision2/interfaces.mp4)

~~~pythonplot Evolution of the kinetic energy
import numpy as np
import matplotlib.pyplot as plt
	
plt.figure()
data1 = np.loadtxt('out', delimiter=' ', usecols=[0,1])
plt.plot(data1[:,0], data1[:,1], ls='--', lw=1.4, label='three-phase w/o conserving');
plt.legend(loc='best');
plt.xlabel('Time')
plt.ylabel('Kinetic Energy')  
plt.xlim([0,3])
plt.tight_layout()
plt.savefig('plot_KE_vs_time.svg')
~~~

We overload *rho* and *mu* so the properties are computed as 
$$ \hat{\rho} (f_1, f_2) = f_1\,\rho_1 + f_2\,\rho_2 + f_3\,\rho_3, \quad \text{where} \quad f_3=1-f_1-f_2 $$
$$ \hat{\mu} (f_1, f_2) = f_1\,\mu_1 + f_2\,\mu_2 + f_3\,\mu_3, \quad \text{where} \quad f_3=1-f_1-f_2 $$
The reason we choose this equation because the sum of coefficients of $A_i$ is one.

*/

#include "navier-stokes/centered.h"
#ifndef rho
#define rho(f1, f2) (clamp(f1, 0., 1.) * rho1 + clamp(f2, 0., 1.) * rho2 + clamp((1.-f1-f2), 0., 1.) * rho3)
#endif
#ifndef mu
#define mu(f1, f2)  (clamp(f1, 0., 1.) * mu1  + clamp(f2, 0., 1.) * mu2  + clamp((1.-f1-f2), 0., 1.) * mu3)
#endif
#include "three-phase.h"               
#include "tension_three-phase.h"


/** Define mesh and simulation time */
#define EPS 1e-10
#define MINLEVEL 4
#define MAXLEVEL 7
#define TMAX 4.0
#define TMOVIES 0.05

/** Define physical parameters for the problem */
#define REYNOLDS 1460
#define WEBER 0.1
#define DENSITY_RATIO 833
#define VISCOSITY_RATIO 55
#define DROP_RADIUS 1.
#define DROP_VELOCITY 1. 
#define DOMAIN_SIZE 8.0 
#define DENSITY_LIQUID 1.
#define DENSITY_GAS (DENSITY_LIQUID/DENSITY_RATIO)
#define VISCOSITY_LIQUID (DENSITY_LIQUID*DROP_VELOCITY*DROP_RADIUS/REYNOLDS)
#define VISCOSITY_GAS (VISCOSITY_LIQUID/VISCOSITY_RATIO)
#define SURFACE_TENSION 1./WEBER

int main(int argc, char** argv)
{  
  L0 = DOMAIN_SIZE [0];
  DT = 1e-3 [0];
  init_grid(1<<(MAXLEVEL));
  origin(-L0/2.,-L0/2.);
  rho1      = DENSITY_LIQUID;
  rho2      = DENSITY_LIQUID;
  rho3      = DENSITY_GAS;
  mu1       = VISCOSITY_LIQUID;
  mu2       = VISCOSITY_LIQUID;
  mu3       = VISCOSITY_GAS;
  sigma12   = SURFACE_TENSION;
  sigma23   = SURFACE_TENSION;
  sigma13   = SURFACE_TENSION;
  
  foreach_dimension()
    periodic(right);
  
  TOLERANCE = 1e-6 [*];
  run();
}

/** Initializes the domain */
#define CIRCLE(x0,y0,R0) (sq(R0) - sq(x-x0) - sq(y-y0))
event init (i=0){
  int ndroplets = 2;
  double drop_radii[ndroplets];
  coord drop_centers[ndroplets];

  for (int j = 0; j < ndroplets; j++)
    drop_radii[j] = DROP_RADIUS; 

  drop_centers[0] = (coord) { 1.0, 1.0};
  drop_centers[1] = (coord) {-1.0,-1.0};
  
  vertex scalar phi[];
  foreach_vertex(){
    phi[] = CIRCLE( drop_centers[0].x, drop_centers[0].y, drop_radii[0] ) ;
  }
  fractions (phi, f1);

  foreach_vertex(){
    phi[] = CIRCLE( drop_centers[1].x, drop_centers[1].y, drop_radii[1] ) ;
  }
  fractions (phi, f2);

  foreach() {
    u.x[] = -DROP_VELOCITY*sqrt(0.5)*f1[] + DROP_VELOCITY*sqrt(0.5)*f2[];
    u.y[] = -DROP_VELOCITY*sqrt(0.5)*f1[] + DROP_VELOCITY*sqrt(0.5)*f2[];
  }
}
#undef CIRCLE

#if TREE
event adapt(i++){
  adapt_wavelet ((scalar *){p, rhov, u}, (double[]){1e-2, 1e-2, 1e-2, 1e-2}, MAXLEVEL, MINLEVEL);
}
#endif

event logfile(i++; t <= TMAX){
  fprintf(stderr, " residuals: %.5f \t %g %d %d \t %g %d %d \n", t, mgp.resa, mgp.i, mgp.nrelax, mgu.resa, mgu.i, mgu.nrelax);
}

/** Outputs a video of the vorticity with interface fragments.
*/

#include "view.h"
event movie (t = 0,t += TMOVIES) {
  scalar omega[];
  vorticity (u, omega);
  squares(color = "omega");
  for (scalar s in interfaces){
    draw_vof(s.name, map = cool_warm, min = 0, max = 3, lw = 2.);
  }
  save ("vorticity.mp4");
  clear();

  int index = 0;
  float * colors[] = {(float[]){1,0,0},(float[]){0,1,0},(float[]){0,0,1},(float[]){1,1,1}};
  for (scalar s in interfaces)
    draw_vof (s.name, fc = colors[index++], filled = 1, min = 0, max = 3, lw = 2.);
  cells(lw=0.5);
  save ("interfaces.mp4");
  clear();
}

/** Track the evolution of kinetic energy
*/

event compute_energies (i++) {
  double vol = 0., rke = 0., se = 0.;
  foreach (reduction(+:vol) reduction(+:rke) reduction(+:se)) {
    vol += dv();
    foreach_dimension() 
      rke += rho(f1[], f2[])*sq(u.x[])*dv();

    for (scalar s in interfaces)
    if (s[] > EPS && s[] < 1 - EPS) {
      coord  normal = mycs(point, s), parea;
      double alpha = plane_alpha(s[], normal);
      double dS = pow(Delta, dimension - 1) * plane_area_center(normal, alpha, &parea);
      se += sigma12*dS; // only because I used sigma12=sigma23=sigma=13
    }
  }
  rke /= 2.*vol;
  fprintf (stdout,"%g %g %g\n",t, rke, se);
  fflush (stdout);
}

