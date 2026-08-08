/**
# Kelvin-Helmholtz simulation using the streamfunction-vorticity Navier-Stokes solver, and the SFD method for stabilisation

The flow is defined with the velocity $U(y)$ as :
$$ U(y) = U_m + \frac{\Delta U}{2} \tanh{\left( \frac{y}{e} \right)} $$
With $U_m = \frac{U_{top} + U_{bottom}}{2}$

Which gives : 
$$
\left\{
\begin{aligned}
    &\omega = - \frac{\Delta U}{2e} \mathrm{sech}\left( \frac{y}{e} \right)^2\\
    &\psi = -U_m y - \frac{e \Delta U}{2}\log\left(\cosh\left(\frac{y}{e}\right)\right)
\end{aligned}
\right.
$$

The boundary conditions are periodic on the sides.

# Code
*/


#include "navier-stokes/stream.h"
#include "SFD_stream.h"

double U_top    = 10. [1,-1];
double U_bottom = 8. [1,-1];

double e      = 0.01 [1];     // thickness of the initial noise
double k      = 1e-2 [0,-1];  // scale of the initial noise
float  lambda = 3.141592653389793 [-1]; // wave-number of the initial perturbation

float Tend = 4;

/** The SFD is activated at $t = 2$. */
double freq_SFD = 4;
bool SFD_toggle;
event adapt_toggle (i++)
  SFD_toggle = (t >= 2);

int main() {

  L0 = 10. [1];
  origin(0, -L0/2);
  N = 256;

  periodic(left);
 
  run();
}

psi[top]    = dirichlet(-(U_top+U_bottom)/2.*y - e*fabs(U_top-U_bottom)/2.*log(cosh(y/e)));
psi[bottom] = dirichlet(-(U_top+U_bottom)/2.*y - e*fabs(U_top-U_bottom)/2.*log(cosh(y/e)));

event init (i = 0) {

/** The grid is refined where the instability develops */
  refine (fabs(y) < 0.7 && level < 10);

/** Initialisation of the vorticity + small perturbations */
  foreach() {
    omega[] = - fabs(U_top-U_bottom)/(2.*e) * sq(1./cosh(y/e));
    omega[] +=  fabs(y) < e ? k*cos(lambda*x) : 0. [0, -1] ;
  }
  poisson(psi, omega);
}

/** Ouputs */

event output (i += 4; t <= Tend)
  output_ppm (omega, box = {{0,-2.5},{L0,2.5}}, min = -15, max = 0.5, file = "omega.mp4");


event logfile (i++)
  fprintf (stderr, "%d %g %g %g\n", i, t, interpolate(omega, L0/2, 0), interpolate(omega, L0/2, L0/4));

/**
# Visualisations

The SFD is activated at t = 2.

![Animation of the vorticity field](KelvinHelmholtzSFD/omega.mp4)(loop)
*/