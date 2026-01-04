/**
# Impinging Planar Jets — Bertsch *et al.* (2020) & Bongarzone *et al.* (2021)

This simulation reproduces the **planar impinging jets configuration** studied in  
- *Bertsch et al., "Feedback-free microfluidic oscillator with impinging jets", Phys. Rev. Fluids, 2020*  
- *Bongarzone et al., "Impinging planar jets: hysteretic behaviour and origin of the self-sustained oscillations", J. Fluid Mech., 2021*  

In this configuration, two **opposing laminar jets** impinge at the mid-plane of a rectangular domain, forming a stagnation region and recirculation cells.  
At moderate Reynolds numbers, the steady symmetric solution can lose stability and lead to periodic oscillations even **without any external feedback** — a purely hydrodynamic oscillator.

This case studies that geometry and basic setup in a 2D viscous, incompressible flow using the *Basilisk* framework. **Many thanks to Timothée Salamon for the collaboration and the linear stability analysis.**

## Simulation Overview

We solve the incompressible Navier–Stokes equations with a passive tracer field `s` to visualize the jet interaction and recirculating structures.

Parameters:
- $Re = 25$ (based on jet radius $R_d$ and inlet velocity $U_0$)
- $R_d = 0.03$
- $U_0 = 1$
- Domain size $L_0 = 1$

Slip/outflow conditions are imposed on the left and right sides; two opposing inflows are prescribed at the bottom and top boundaries.

*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"

scalar s[];
scalar * tracers = {s};

double U0;
double R_d;
double Re;

face vector muv[];
FILE * fpmax;

int main() {

  Re = 25;         // moderate Reynolds number, laminar regime
  R_d = 0.03;      // half jet width
  L0  = 1.0;       // domain size
  U0  = 1.0;       // reference velocity

  TOLERANCE = 1e-3 [*];

  
  //- **Bottom boundary:** upward jet (positive y-direction)
  //- **Top boundary:** downward jet (negative y-direction)
  //- **Left & Right:** open / outflow conditions
  
  //Each jet is defined by a rectangular inlet of half-width `R_d`.
  

  u.n[bottom] = dirichlet ( U0 * (x > -R_d && x < R_d) );
  u.t[bottom] = dirichlet (0.);
  s[bottom]   = dirichlet ( 2 * (x > -R_d && x < R_d) );

  u.n[top] = dirichlet ( -U0 * (x > -R_d && x < R_d) );
  u.t[top] = dirichlet (0.);
  s[top]   = dirichlet ( -2 * (x > -R_d && x < R_d) );

  u.n[left] = neumann(0);
  p[left]   = dirichlet(0);

  u.n[right] = neumann(0);
  p[right]   = dirichlet(0);

  origin(-L0/2, 0);
  N = 128;
  init_grid(N);
  
  char param_dim[80];
  sprintf(param_dim, "param_dim.txt");
  FILE * fparam = fopen(param_dim, "w");
  fprintf(fparam, "%g %g %g %d\n", L0, R_d, U0, N);
  fclose(fparam);

  mu = muv;
  fpmax = fopen("log.dat", "w"); 

  run();
}

/**
## Initialization

We initialize the tracer `s` to distinguish the top and bottom flow regions visually.
*/

event init (t = 0) {
  foreach()
    s[] = (y < L0/2);  // simple step distribution
}

/**
## Viscosity Field

We enforce a uniform kinematic viscosity corresponding to the chosen Reynolds number.
*/

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[] * U0 * R_d / Re;
}

event logfile (i++) { 
  fprintf(stderr, "%d %g\n", i, t); 
  fprintf(fpmax, "%d %g\n", i, t);
}

/**
## Visualization Outputs

We output the vertical velocity component and the tracer field for visualization.
*/

event ppm_output (t = 0; t += 0.1; t <= 100) {

  char name1[80];
  sprintf(name1, "uY.mp4");
  output_ppm(u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

  char name2[80];
  sprintf(name2, "s.mp4");
  output_ppm(s, file = name2, n = 512, min = -2, max = 2, linear = true);
}

/**
![Passive tracer field](imp_jets/s.mp4)
![Vertical velocity field](imp_jets/uY.mp4)

We observe the interaction region where the opposing jets collide at mid-height,
forming recirculating cells and a stagnation point.
For Reynolds numbers near the onset threshold ($Re \approx 25$),
small asymmetries can grow, leading to a **self-sustained oscillation** of the impingement point,
as reported by Bongarzone *et al.* (2021) and Bertsch *et al.* (2020).
*/

event profile (t = end) {
  printf("-----END-----\n");
}
/**
## Results and interpretations

A clear unstable periodic mode is taking place in the channel. We can perform, as done in the articles, a linear stability study. Here are some quick results obtained by Timothée Salamon.

<img src="https://lh3.googleusercontent.com/sitesv/AAzXCkcqvbCUkOffg6Cb71PZ-FUt3itkyoR65ezSEBHJx8KpFiNY4PF4WtVgeY78usICMe1MXEquUjNsfqzHU7gODGDoigCvp42hMA0HmQOPO59MhtXs6ENcAQhjBy7RFavnzqmGbvuvez1h-e4x5CTOJrfXX-1U6zhwX1NR84PUZA6xOGhQiVzBx7vm2jcOwKV0smzgkXk7HbA2bjMp2tgaVbs0eXUmO6QkjTPnsH0=w1280" alt="drawing" width="400" style="display:block; margin-left:auto; margin-right:auto;"/>

<p align="center">*Mode 1 obtained by linear stability analysis*</p>

Mode 1 is the unstable mode that we are mostly seeing. It has a exponential growth rate of $\approx 2 \times 10^{-2}$ s$^{-1}$ and a pulsation of $\approx 1.5 \times 10^{-2}$ s$^{-1}$ which is the one we see on the video.

<img src="https://lh3.googleusercontent.com/sitesv/AAzXCkdi6eoBX5gGY3ngN4uupuUyWtF804ik3P5mAcVn4r0xouwqYt1HMFEYSy028cDZE3SBxQ7cHhk0OdU9OV3MEzyJtjmzytyjIawiNITgIG7KRSFwQtjrQ4MO9xHGD_k8VbVdpBRl_T4OL89AdIQkWQ9kpSrW3ImsGMuvR7cW2Pt-XaftaHVa4pmiV4p-R0OFKLq949V4qC3lIYcApgbSR2IyJlip838luccp=w1280" alt="drawing" width="400" style="display:block; margin-left:auto; margin-right:auto;"/>

<p align="center">*Mode 2 obtained by linear stability analysis*</p>

<img src="https://lh3.googleusercontent.com/sitesv/AAzXCkezU6YfWGV5TUcuIntSCjpOO6kfpUd3n6hrvR4JaWN7y_5BxHcM3Gm8L84mPgp1aR6ksXN8121zaJ7clG8hswJg6NAOiO_qsRiPXU-MOJpyY7gufpebZgs4tih_m4wYWC5ERFS9FsEZOGfSuLYX6jY62uql4DkTbFlaMOQ--BidcbH-zEAmtOv-uq8KxQHBD-iHo1dLIRxzUgDssg6W3_h-ttJfAYEzlUgS98E=w1280" alt="drawing" width="400" style="display:block; margin-left:auto; margin-right:auto;"/>

<p align="center">*Mode 3 obtained by linear stability analysis*</p>

Mode 2 and Mode 3 are stationnary mode (meaning the pulsation $\omega$ = 0). They are very close one to another in the eigenvalue plot (they are circled), so maybe they are in reality the same mode and it is the mesh discretization that creates two separated spurious modes.

<img src="https://lh3.googleusercontent.com/sitesv/AAzXCkcW6H7OylpTrYbU8L3Ytvw9GUnhGFJWGdCRy4Bjwx_wOqJDTOPEu_FvJI_qNkgPnMxp2BYJy-9IQk_OIYw4TA6gwhVm5pQVmlbuMC7Hl8ux-LW5QPPvacivTTVvgt8ntvgzgyzY0TuSZFJIKBTid2DRvcWC02gVrGon-FyJ1ySOmJ9loMR-ZPaypGH-bleibnxkC0xvwTOpFdaX1QpEeJx5HLZe3NLKkAAE=w1280" alt="drawing" width="400" style="display:block; margin-left:auto; margin-right:auto;"/>

<p align="center">*Eigenvalue map*</p>

*/
