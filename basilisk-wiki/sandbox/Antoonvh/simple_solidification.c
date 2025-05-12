/**
![Ice also has its downside, here it may exchange heat with
water. Image via [Wageningen Univeristy and
Research](https://www.wur.nl/nl/show/Penguins-give-underwater-show-in-front-of-ice-net-in-Antarctic.htm)](https://www.wur.nl/upload_mm/4/4/b/8d7514c0-dfbb-4e53-b963-4a13629b5477_Franeker-snapshot001.jpg)

# Simple solidification

This page concerns a simple model for solidification and melt.

## Setup

Consider a fluid with that freezes at $T = 0$ with an initial
temperature gradient $T = ay$, with $a > 0$ being a constant
temperature gradient and $y$ the vertical coordinate. A warm dipolar
vortex of size $R$ is targeted at the fluid-solid interface with
downward propagation velocity $U$. It is initialized at a distance
$y_o$ away from the fluid-solid interface. At the fluid-solid
interface, a boundary layer may emerge that owes its existense to the
fluid's viscosity ($\nu$). Futhermore, due to the transport of heat,
the solid may melt or solidify. There is diffusion of heat via the
medium's constant diffusivity ($\kappa$). One may identify three
dimensionless groups:

$$Re = \frac{UR}{\nu}, $$

$$\Pi = \frac{y_o}{R}, $$

$$Pr= \frac{\kappa}{\nu}, $$

and I *expect* that the exact value of $a$ is not relevant for the
dynamics. We choose $Re = 1000$, $\Pi = 4$ and study the effects of
various values for $Pr$. Unitary values are chosen for $R, U, a$.
 */

#include "navier-stokes/centered.h"
#include "tracer.h"

scalar T[], * tracers = {T};
T[top]    = dirichlet (y);
T[bottom] = dirichlet (y);

double Re = 1000.;
double Pi = 4.;
double Pr;
int maxlevel = 10;

event solid (i++) {
  foreach() {
    if (T[] < 0)
      foreach_dimension()
	u.x[] = 0;
  }
}

int main () {
  L0 = 20;
  X0 = Y0 = -L0/2.;
  N = 256;
  for (Pr = 1; Pr <= 100; Pr *= 10)
    run();
}

/**
## Initialization

The initial Lamb-Chaplygin vortex is placed in the warm part of the
fluid.
*/
double k = 3.83170597;
double xo = 0.1;
double yo = 5;
#define RAD (pow(pow((x - xo), 2)+(pow((y - yo), 2)), 0.5))
#define ST (-(x - xo)/RAD)
event init (t = 0) {
  const face vector muc[] = {1./Re, 1./Re};
  mu = muc;
  refine (RAD < 2.0 && level <= 9);
  refine (RAD < 1.0 && level <= 10);
  scalar psi[]; //Stream function
  foreach() 
    psi[] = ((RAD > 1)*((1/RAD))*ST) + ((RAD < 1)*((-2*j1(k*RAD)*ST/(k*j0(k))) + (RAD*ST)));
  boundary ({psi});
  foreach() {
    u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta));
    u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta);
    T[] = y; 
  }
}

#include "diffusion.h"
event tracer_diffusion (i++) {
  const face vector kappa[] = {Pr/Re, Pr/Re};
  diffusion (T, dt, kappa);
}

event adapt (i++) 
  adapt_wavelet ({T, u}, (double[]){0.1, 0.025, 0.025}, maxlevel);
/**
## Output

For each run, movies are generated revealing the evolution of the
temperature and the vorticity fields. This is supplemented with a
movie of the adaptive mesh-structure.

![Vorticity, $Pr = 1, 10, 100$](simple_solidification/vor.mp4)

![Temperature $T$, $Pr = 1, 10, 100$](simple_solidification/T.mp4)

![Level of refinement, $Pr = 1, 10, 100$](simple_solidification/level.mp4)
*/
event movie (t += 0.1) {
  scalar m[], omega[], lev[];
  foreach() {
    m[] = (T[] > 0) - 0.5;
    lev[] = level;
  }
  vorticity (u, omega);
  
  output_ppm (omega, mask = m, file = "vor.mp4", n = 400);
  output_ppm (T, file = "T.mp4", n = 400);
  output_ppm (lev, file = "level.mp4", min = 2, max = maxlevel, n = 400);
}

event stop (t = 15);
