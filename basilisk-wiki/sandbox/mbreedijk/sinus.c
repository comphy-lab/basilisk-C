/**
# Convection in a stably stratified atmosphere

We consider a stably stratified atmosphere with a constant Brunt
Vaisiala frequency $N$. We use buoyancy as our thermodynamic variable.
The initial vertical profile is,

$$b(y, t = 0) = N^2y.$$

At the surface ($y = 0$), the mean buoyancy is equal to $b(0, t = 0)$
but with a sinusoidal dependency in space,

$$b(x, y = 0) = b_0 \mathrm{sin}(k x),$$

where $k$ is the wave number $[L^{-1}]$ and $b_0$ a surface-buoyancy
amplitude scale. We can construct a relevant lengthscale for the
convective motions ($\mathcal{L}$).

$$\mathcal{L} = \frac{b_0}{N^2}.$$

This gives rise to a dimensionless group, comparing the
periodicity-length scale ($\propto k^{-1}$) and $\mathcal{L}$,

$$ \Pi_1 = \mathcal{L}k.$$

In order to make the simulations feasible and create artificial
dissipation of second order moments, we consider an atmosphere with
hyper viscousity $\nu$ and diffusivity $\kappa$. Intoducing, 

$$\Pi_2 = \frac{b_0^2 N^3}{\nu},$$

$$\Pi_3 = \frac{\nu}{\kappa}.$$
 
For the set of 5 parameters ($b_0, N, k, \nu, \kappa$) with 2
base-units ($L, T$), $\Pi_1, \Pi_2$ and $\Pi_3$ form a ortogonal group
of dimensionless numbers that fully describe the system. We choose

$$\Pi_1 = 3,$$
$$\Pi_2 = 2000,$$
$$\Pi_3 = 1.$$

## Setup 

We want to solve the two-dimensional Navier-Stokes equations for
incompressible fluid flow with a buoyancy tracer.
 */
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
scalar b[], * tracers = {b};
/**
We set unitary values for $b_0$ and $N$, and declare the other
parameters.
 */
double b0 = 1, Nbv = 1, k, nuv, kappav;
double Pi1 = 3, Pi2 = 2000, Pi3 = 1;
/**
The boundary conditions are of the Dirichlet and Neumann type.
 */
b[bottom] = dirichlet (b0*sin(k*x));
b[top]    = neumann (sq(Nbv));
u.t[bottom] = dirichlet (0); //no slip
/**
The domain size is chosen $L_0 = \frac{4 \pi}{k}$ so that we can use
periodic boundary conditions in the horizontal (`left/right`)
direction. Further, we must link the gravity-acceleration vector.
 */
face vector av[]; 
int main() {
  periodic (left);
  k = Pi1/(b0/sq(Nbv));
  nuv = (b0*cube(Nbv))/Pi2;
  kappav = nuv;
  L0 = 4*pi/k;
  const face vector muc[] = {nuv, nuv};
  // Link fields
  mu = muc;
  a  = av;
  run();
}
/**
The stratification is initialized, tracer-diffusion is switched on,
the gravity vector is computed and the grid is adapted in the
following events.
 */
event init (t = 0) {
  DT = 0.1/Nbv;
  foreach() 
    b[] = sq(Nbv)*y;
  boundary ({b});
}

event tracer_diffusion (i++) {
  const face vector kappa[] = {kappav, kappav};
  diffusion (b, dt, kappa);
}

event acceleration (i++) {
  coord g = {0, 1};
  foreach_face()
    av.x[] = g.x*face_value(b, 0); 
}

event adapt (i++) {
  adapt_wavelet ({u, b}, (double[]){0.01, 0.01, b0/50.}, 8);
}
/**
### Damping

At the top of the domain we should damp internal waves.
*/
 
event damp (i++) {
  foreach() {
    if (y > L0/2) {
      foreach_dimension()
	u.x[] *= exp(-(y - L0/2)/5.);
    }
  }
  boundary ((scalar*){u});
}

/**
### Output

A few movies
 */
event mov (t += 1) {
  output_ppm (b, file = "b.mp4", n = 300,
	      min = -1, max = 1.5, linear = true);
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (lev, file = "l.mp4", n = 300, min = 1, max = 8.5);
}

event stop (t = 500);

/**
### Results

![Buoyancy](sinus/b.mp4)

![Grid](sinus/l.mp4)
 */ 
