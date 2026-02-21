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
