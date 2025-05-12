/**
# Description
<figure>
<p align="center">
  <img src="https://www.dropbox.com/s/4s20uxvgfbxj4xz/InitialConditionDropOnDropImpact.png?dl=1" width="25%">
  <figcaption><p align="center">We want to reproduce the shape of this drop. The orange contour is from the numerical simulation, and the image is from experiments [here](http://basilisk.fr/sandbox/vatsal/DropOnDropImpact/dropOnDropImpact.c). Click [here](https://www.dropbox.com/s/4s20uxvgfbxj4xz/InitialConditionDropOnDropImpact.png?dl=0) if the above figure is not displayed properly on your browser.</figcaption>
</figure>
For Bond numbers $> 0.1$, assuming the sessile drop to be spherical is inaccurate. One can also see this in the experimental image above. Another example is [(here)](http://basilisk.fr/sandbox/Antoonvh/rest.c). I wanted to be as close to the experiments as possible. Hence, the code

# Numerical Code
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
/**
We use a modified adapt-wavelet algorithm available [(here)](http://basilisk.fr/sandbox/pairetti/bag_mode/adapt_wavelet_limited.h). It is written by *César Pairetti* (Thanks :)). We use to ensure that refinement is higher near the substrate.
*/
#include "pairetti/bag_mode/adapt_wavelet_limited.h"

int MAXlevel = 6; // maximum level: This is increased in time.
#define MINlevel 4 // minimum level
#define Ldomain 4
#define tmax 25
#define tsnap (0.05)
// Error tolerances
#define fErr (1e-3) // error tolerance in VOF
#define KErr (1e-4) // error tolerance in KAPPA
#define OmegaErr (1e-3) // error tolerances in vorticity

#define Mu21 (6e-3) // viscosity ratio between the gas and the liquid
#define Rho21 (1./770) // density ratio between the gas and the liquid
char comm[80], FacetName[80], nameOut[80];

/** Boundary Conditions
Substrate is superamphiphobic and has the no-slip condition for velocity.
*/
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

int main(){
  /**
  Navier Stokes equation for this case:
  $$
  \partial_tU_i+\nabla\cdot(U_iU_j) =
  \frac{1}{\hat{\rho}}\left(-\nabla( p - Bo g_jX_j) + Oh\nabla\cdot(2\hat{\mu}D_{ij}) + (\kappa - Bo\Delta\rho g_jX_j)\delta_sn_i\right) + Bog_i
  $$
  The $\hat{\rho}$ and $\hat{\mu}$ are the VoF equations to calculate properties, given by:
  $$
  \hat{A} = (f_1+f_2) + (1-f_1-f_2)\frac{A_g}{A_l}
  $$
  */
  L0=Ldomain;
  X0=0.; Y0=0.;
  init_grid (1 << (4));
  /**
  Bond number $Bo$: measure between Gravity and surface tension.
  $$ Bo = \frac{\rho_lgR^2}{\gamma} $$
  We use [reduced.h](http://basilisk.fr/src/reduced.h) implementation for gravity.
  */
  double Bond = 0.308;
  G.x = -Bond;
  /**
  Velocity scale as the intertial-capillary velocity,
  $$ U_\gamma = \sqrt{\frac{\gamma}{\rho_l R}} $$
  */
  f.sigma = 1.;
  /**
  Ohnesorge number $Oh$: measure between surface tension and viscous forces.
  $$ Oh = \frac{\mu_l}{\sqrt{\rho_l\gamma R}} = 1 $$
  Oh number is assumed to be 1 for faster convergence.
  */
  mu1 = 1.; mu2 = Mu21;
  rho1 = 1.; rho2 = Rho21;
  run();
}

event init (t=0){
  sprintf (comm, "mkdir -p intermediate%5.4f",  fabs(G.x));
  system(comm);
  refine (sq(x-1.015625) + sq(y) < sq(1.04) && level < MAXlevel+2);
  fraction (f, 1. - sq(x-1.015625) - sq(y));
}

/**
Write log and check for convergence
*/
event logWriting (i=i+25) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += 2*pi*y*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke MaxLevel\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g %d\n", i, dt, t, ke, MAXlevel);
  /**
  We stop the simulation when Kinetic Energy (overall) is below $10^{-6}$.
  */
  if(i>1e4 && ke<1e-6){
    fprintf(ferr, "Convergence Reached for Bo = %f\n",-G.x);
    dump (file = "dumpConverged");
    return 1;
  }
}

/**
This event is specific to César's adapt_wavelet_limited. Near the substrate, we refine the grid one level higher than the rest of the domain.
*/
int refRegion(double x, double y, double z){
  return (x < 0.128 ? MAXlevel+2 : x < 0.256 ? MAXlevel+1 : MAXlevel);
}

/**
## Adaptive Mesh Refinement
*/
scalar KAPPA[], omega[];
event adapt(i++){
  MAXlevel = t < 5.0 ? 6 : t < 10.0 ? 7 : 8;
  curvature(f, KAPPA);
  vorticity (u, omega);
  foreach(){
    omega[] *= (f[]);
  }
  boundary((scalar *){KAPPA, omega});
  adapt_wavelet_limited ((scalar *){f, KAPPA, omega},
     (double[]){fErr, KErr, OmegaErr},
     refRegion, MINlevel);
}

/**
## Dumping snapshots
*/
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = "dump");
  sprintf (nameOut, "intermediate%5.4f/snapshot-%5.4f", fabs(G.x), t);
  dump (file = nameOut);
}


/**
## Running the code

Use the following procedure:

~~~bash
#!/bin/bash
mkdir intermediate
qcc -fopenmp -O2 -Wall DropDeposition.c -o DropDeposition -lm
export OMP_NUM_THREADS=8
./DropDeposition
~~~

# Output and Results
The post-processing codes and simulation data are available at: [PostProcess](https://www.dropbox.com/sh/7thev8gv3k7f87b/AAC1ei8rVcFVZ2lQY0HGto7Ba?dl=0)

## The process
<p align="center">
<video width="50%" controls>
  <source src="https://www.dropbox.com/s/ztoiwu49qq99v0q/InitialCondition.mp4?dl=1" type="video/mp4">
  <caption><p align="center">Video available [here](https://www.dropbox.com/s/ztoiwu49qq99v0q/InitialCondition.mp4?dl=0) as well.</caption>
</video>

# Usuage
## Example
* [Drop on drop impact case](http://basilisk.fr/sandbox/vatsal/DropOnDropImpact/dropOnDropImpact.c)
*/
