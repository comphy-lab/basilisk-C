/**
# Harmonic decomposition in a periodic oscillatory flow

This example illustrates the use of [harmonic.h](/src/harmonic.h) on a simple, periodic, 2D oscillatory flow.

We impose a time-periodic inlet velocity on both the left and right sides of a 2D domain.  
The goal is to demonstrate how the harmonic decomposition separates different frequency components of the velocity field.

## Simulation
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "harmonic.h"

double Re;
double U0;

double omegas[3];

const double omega0 = 2.*pi;

face vector muv[];

FILE * fpmax;

int main() {

  Re = 1.;
  U0 = 0.1;
  L0 = 1.; 

  TOLERANCE = 1e-3 [*]; 

  /**
  Periodic-like oscillations are imposed on both left and right sides.
  */

  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);

  u.n[left]  = dirichlet(U0*sin(omega0*1.0*t) * 4.*y/L0 * (1. - y/L0));
  u.n[right] = dirichlet(U0*sin(omega0*1.0*t) * 4.*y/L0 * (1. - y/L0));
  u.t[right] = neumann(0.);

  N = 128;
  origin (0, 0);  // set the origin
  init_grid(N);

  fpmax = fopen("log.dat", "w");

  mu = muv;

  run();
}

/**
We define the two frequencies that will be analysed using `harmonic.h`.
*/
event init (t = 0) {
  omegas[0] = omega0*1.0; 
  omegas[1] = omega0*1.5;  
  omegas[2] = 0.0;           // sentinel value for harmonic.h
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[] * U0 * L0 / Re;
}

event logfile (i++) { 
  fprintf (stderr, "%d %g\n", i, t);
  fprintf (fpmax, "%d %g\n", i, t);
}

/**
The function `harmonic_decomposition()` (from [harmonic.h](/src/harmonic.h)) performs a least-squares regression of a scalar or vector field — here, $u_x$ — over a set of prescribed angular frequencies $\omega_i$.

It fits the temporal signal at each grid point as:
$$
\tilde{u_x}(x,y,t) = Z(x,y) + \sum_i [A_i(x,y)\cos(\omega_i t) + B_i(x,y)\sin(\omega_i t)]
$$

where $A_i$ and $B_i$ are spatial fields representing the local amplitude and phase of each mode.
*/

event harmonic_analysis (t += 0.01; t <= 5) {
  harmonic_decomposition(u.x, t, omegas);
}

event profile (t = end) {
  printf ("----- END -----\n");
}

event ppm_output (t = 0; t += 0.01; t <= 5) {

  char name[80];
  sprintf (name, "uX.mp4");
  output_ppm (u.x, file = name, n = 512, min = -0.1, max = 0.1, linear = true);

  // amplitude for mode 1 (ω₁)
  scalar amp_mode_1[];
  int k1 = 0;
  foreach()
    amp_mode_1[] = sqrt(sq(val(u.x.harmonic.A[k1])) + sq(val(u.x.harmonic.B[k1])));
  boundary ({amp_mode_1});
  char name_u_1[80];
  sprintf(name_u_1, "u_amp_1.mp4");
  output_ppm(amp_mode_1, file = name_u_1, n = 512, min = 0, max = 0.1, linear = true);

  // amplitude for mode 2 (ω₂)
  scalar amp_mode_2[];
  int k2 = 1;
  foreach()
    amp_mode_2[] = sqrt(sq(val(u.x.harmonic.A[k2])) + sq(val(u.x.harmonic.B[k2])));
  boundary ({amp_mode_2});
  char name_u_2[80];
  sprintf(name_u_2, "u_amp_2.mp4");
  output_ppm(amp_mode_2, file = name_u_2, n = 512, min = 0, max = 0.1, linear = true);
}

/**
## Visualization

![Horizontal velocity field](test_harmonic/uX.mp4)
![Projection of the velocity on mode #1](test_harmonic/u_amp_1.mp4)
![Projection of the velocity on mode #2](test_harmonic/u_amp_2.mp4)

*/

/**
# Interpretation of the harmonic decomposition

The *harmonic analysis* produces maps of the amplitude associated with each tested frequency $\omega_i$.

A strong, nearly time-independent amplitude field indicates that the corresponding frequency is **dominant** in the flow — meaning the local velocity oscillates periodically with that frequency.

For instance, in this example:

- The mode at $\omega_1 = 2\pi$ corresponds to the imposed oscillation from the boundaries. Its amplitude field shows where the oscillation penetrates into the domain.
- The mode at $\omega_2 = 3\pi$ corresponds to a higher harmonic. Its amplitude rapidly goes to nearly zero.
*/
