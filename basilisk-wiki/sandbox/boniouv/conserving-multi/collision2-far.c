/**
# Collision of 2 droplets in a periodic box

This test case aims to validate the extension of 
[the momentum conservation algorithm](/src/navier-stokes/centered.h) 
to the [multi-VOF method](no-coalescence.h). 

This case is made non-dimensional using the drop radius $R$ and drop velocity
$U$. The problem parameters are
$$
\mathrm{Re} = 1460, \quad 
\mathrm{We} = 0.1, \quad 
\rho^* = 833, \quad 
\mu^* = 55 
$$

At $t=0$, the droplets are far enough so [no-coalescence.h](no-coalescence.h)
is not active. The combination of 
```
#include "conserving-multi.h" // Remove this to go back to standard scheme
#include "no-coalescence.h"
``` 
and 
```
#include "navier-stokes/conserving.h"
#include "no-coalescence.h"          
```
with these initial conditions results in `f[]` being reset.  

![Results from `collision2-far.c`](collision2-far/interfaces.mp4)


*/


#include "navier-stokes/centered.h"
#include "two-phase-multi.h"
#include "conserving-multi.h"
#include "tension.h"
#include "no-coalescence.h"
       
#include "view.h"
#include "../misc/fractions_init.h"

/** Define mesh and time */
#define EPS 1e-10
#define MINLEVEL 4
#define MAXLEVEL 7
#define TMAX 3.
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
  rho2      = DENSITY_GAS;
  mu1       = VISCOSITY_LIQUID;
  mu2       = VISCOSITY_GAS;
  f.sigma   = SURFACE_TENSION;
  
  foreach_dimension()
    periodic(right);
  
  TOLERANCE = 1e-6 [*];
  run();
}

/** Initializes the domain */
#define CIRCLE(x0,y0,R0) (sq(x-x0) + sq(y-y0) <= sq(R0))
event init (t=0){
  int ndroplets = 2;
  double drop_radii[ndroplets];
  coord drop_centers[ndroplets];

  for (int j = 0; j < ndroplets; j++)
    drop_radii[j] = DROP_RADIUS; 

  drop_centers[0] = (coord) { 1.0, 1.0};
  drop_centers[1] = (coord) {-1.0,-1.0};

  foreach() {
    u.x[] = 0.;
    u.y[] = 0.;
    
    coord center = {x,y};
    f[] = refine_frac(center, drop_centers, drop_radii, ndroplets, Delta, 1);
    
    if CIRCLE(drop_centers[0].x, drop_centers[0].y, drop_radii[0]){
      u.x[] = -DROP_VELOCITY*sqrt(0.5)*f[];
      u.y[] = -DROP_VELOCITY*sqrt(0.5)*f[];
    } else if CIRCLE(drop_centers[1].x, drop_centers[1].y, drop_radii[1]){
      u.x[] = DROP_VELOCITY*sqrt(0.5)*f[];
      u.y[] = DROP_VELOCITY*sqrt(0.5)*f[];
    }
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

/** Track the evolution of kinetic energy and surface energy.
*/

event compute_energies (i++) {
  double vol = 0., rke = 0., se = 0.;
  foreach (reduction(+:vol) reduction(+:rke) reduction(+:se)) {
    vol += dv();

    foreach_dimension() 
      rke += rho(f[])*sq(u.x[])*dv();

    for (scalar s in interfaces)
      if (s[] > EPS && s[] < 1 - EPS) {
        coord  normal = mycs(point, s), parea;
        double alpha = plane_alpha(s[], normal);
        double dS = pow(Delta, dimension - 1) * plane_area_center(normal, alpha, &parea);
        se += f.sigma*dS;
      }
  }
  rke /= 2.*vol;
  se /= vol;
  fprintf (stdout,"%g %g %g\n",t, rke, se);
  fflush (stdout);
}

