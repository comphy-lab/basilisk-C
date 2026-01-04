/**
# Leidenfrost Droplet - Fixed Flux

The [Leidenfrost effect](https://en.wikipedia.org/wiki/Leidenfrost_effect)
is a physical phenomenon in which a liquid, close to a solid
surface, evaporates producing an insulating vapor layer which
keeps the liquid droplet away from the solid surface.

## Choice of the Navier-Stokes Equations Solver

The numerical simulation of this phenomena is not straight-forward, the main
complication is that the model to obtain the liquid (extended) velocity for the
transport of the volume fraction must consider the influence of the Stefan flow.
Initially, the droplet moves toward the solid surface. The evaporation process
creates a velocity field which points radially outward from the droplet and,
when the drop is sufficiently close to the solid surface, the Stefan flow pushes
the droplet away from the wall.

To achieve this coupling between the transport of the liquid phase and the
Stefan flow we can't use the
[navier-stokes/low-mach.h](/sandbox/ecipriano/src/navier-stokes/low-mach.h)
model with 2 velocity fields, because this model decouples the advection
velocity from the field velocity which includes the Stefan flow. Although this
approach is beneficial when dealing with static droplets with strong density
ratio, it is not adequate for Leidenfrost effect simulations. We use the
[navier-stokes/velocity-jump.h](/sandbox/ecipriano/src/navier-stokes/velocity-jump.h)
approach instead, which limits oscillations in the velocity field, and it
couples the effect of the Stefan flow with the advection of the liquid phase.

## Results

From the video we observe that the droplet never touches the boundary. When the
the drop approaches the wall, we observe the formation of the vapor layer which
pushes the droplet away from the wall.
Depending on the physical properties used, the vapor layer thickness can
become very small, which requires a very fine grid in order to be able to resolve
this effect and to avoid the contact of the droplet with the wall.
If we use a more complex phase change model we need to ensure that the vaporization
rate is not too small, by setting a high value of wall and environment temperature.

![Evolution of the velocity field and the interface position](leidenfrostfixedflux/movie.mp4)
*/

#include "navier-stokes/velocity-jump.h"
#include "two-phase.h"
#include "tension.h"
#include "phasechange.h"
#include "fixedflux.h"
#include "recoil.h"
#include "view.h"

/**
### Boundary Conditions

We set a no-slip boundary condition on the solid wall (left), and outflow
boundary conditions on the opposite side of the domain (right). Symmetry
everywhere else. */

u.n[left] = dirichlet (0.);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

/**
We set the maximum level of refinement, the vaporization rate per unit of
surface, the initial droplet diameter, and the initial droplet velocity.

For the data used in this simulation, the Weber number is:

$$
  We = \dfrac{\rho v_0^2 D}{\sigma} = 4
$$
*/

int maxlevel;
double D0 = 1.e-3;
double v0 = -0.2;

int main (void) {

  /**
  We set the physical properties of the simulation. */

  rho1 = 200., rho2 = 5.;
  mu1 = 1.e-3, mu2 = 1.e-5;
  mEvapVal = -0.8;

  f.sigma = 0.002;

  /**
  We change the length of the domain and we run the
  simulation. */

  L0 = 5.*D0;

  /**
  We use two different velocity fields for the GFM/velocity-jump approach. The
  solution of the temperature and species equations is skipped. */

  nv = 2;
  pcm.isothermal = true;
  pcm.isomassfrac = true;

  for (maxlevel = 7; maxlevel <= 7; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

/**
The droplet is initialized at the coordinate $(x,y) = (2D_0,0)$. */

#define circle(x,y,R) (sq(R) - sq(x - 2.*D0) - sq(y))

event init (i = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  /**
  We initialize the velocity field depending on the model used
  for the solution of the Navier-Stokes equations. */

  foreach()
    for (vector u in ulist)
      u.x[] = v0*f[];
}

/**
We adapt the grid according to the position of the interface and the
velocity field. */

event adapt (i++) {
  adapt_wavelet_leave_interface ({ur.x,ur.y}, {f},
      (double[]){1.e-2,1.e-2}, maxlevel, 1);
}

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

event movie (t += 0.0001, t <= 0.05) {
  clear();
  view (theta=0., phi=0., psi=-pi/2.,
      tx = 0., ty = -0.26, fov = 10.);

  draw_vof ("f", lw = 1.5);
  squares ("ur.x", spread = 0);
  mirror ({0.,1.}) {
    draw_vof ("f", lw = 1.5);
    vectors ("ur", scale = 2.e-5);
  }
  save ("movie.mp4");
}
