/**
# Expansion of a Nucleated Bubble

A spherical bubble is initialized on the lower-left corner of and AXI domanin.
Similarly to the Scriven problem, the bubble is expanded due to the temperature
gradient between the interface of the bubble, which remains at saturation
temperature, and the superheated liquid environment. However, in this case we
consider the left side of the boundary as a solid wall and we impose a contact
angle, which influences the bubble growth process.

This is the simplest nucleate boiling configuration, where no microlayer
modelling and no heat exchange with the solid wall is condered.

![Evolution of the interface, the temperature field, and the grid refinement](bubblecontact/movie.mp4)(width="100%")
*/

/**
## Simulation Setup

We use the centered Navier--Stokes equations solver with volumetric source in
the projection step. The phase change is directly included using the boiling
module, which sets the best (default) configuration for boiling problems. Many
features of the phase change (boiling) model can be modified directly in this
file without changing the source code, using the phase change model object
`pcm`. Compiling with `-DJUMP=1` changes the Navier--Stokes solver to the
velocity-jump formulation, which employs a GFM approach to set the interface
velocity jump. A fixed contact angle is imposed on the left side of the domain.
*/

#include "axi.h"
#if JUMP
# include "navier-stokes/velocity-jump.h"
#else
# include "navier-stokes/low-mach.h"
#endif
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "boiling.h"
#include "view.h"

/**
### Boundary conditions

Outflow boundary conditions for velocity and pressure are imposed on the top and
right walls. On the left wall, no-slip conditions are imposed.*/

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

u.n[left] = dirichlet (0.);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

double Tsat, Tbulk;
T[top] = dirichlet (Tbulk);
T[right] = dirichlet (Tbulk);

/**
To set the contact angle, we allocate a [height-function field](/src/heights.h)
and set the contact angle boundary condition on its tangential component. */

vector h[];
double theta0 = 40.;
h.t[left] = contact_angle (theta0*pi/180.);

/**
### Problem Data

We declare the maximum and minimum levels of refinement, the initial bubble
radius, the growth constant, and additional variable for post-processing. */

int maxlevel = 7, minlevel = 3;
double R0 = 0.1;
double V0, tshift;
double betaGrowth = 0.7659010001953077;

/**
### Initial Temperature

We use gsl to compute the integral required to obtain the Scriven temperature
profile. */

#include <gsl/gsl_integration.h>
#pragma autolink -lgsl -lgslcblas

double intfun (double x, void * params) {
  double beta = *(double *) params;
  return exp(-sq(beta)*(pow(1. - x, -2.) - 2.*(1. - rho2/rho1)*x - 1 ));
}

double tempsol (double r, double R) {
  double result, error;
  double beta = betaGrowth;
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &intfun;
  F.params = &beta;
  gsl_integration_qags (&F, 1.-R/r, 1., 1.e-9, 1.e-5, 1000,
                        w, &result, &error);
  gsl_integration_workspace_free (w);
  return Tbulk - 2.*sq(beta)*(rho2*(dhev + (cp1 - cp2) *
        (Tbulk - Tsat))/rho1/cp1)*result;
}

int main (void) {

  /**
  We set the material properties of the two fluids. In addition to the classic
  Basilisk setup for density and viscosity, we need to define thermal
  properties, such as the thermal conductivity $\lambda$, the heat capacity
  $cp$, and the enthalpy of vaporization $\Delta h_{ev}$. */

  rho1 = 2.5; rho2 = 0.25;
  mu1 = 7.e-3; mu2 = 7.e-4;
  lambda1 = 0.07, lambda2 = 0.007;
  cp1 = 2.5, cp2 = 1.;
  dhev = 100.;

  /**
  The initial bubble temperature and the interface temperature are set to the
  saturation value. */

  Tbulk = 3., Tsat = 1., TIntVal = 1.;
  TL0 = Tbulk, TG0 = TIntVal;

  /**
  We solve two different sets of Navier--Stokes equations, which is fine for
  both the double pressure-velocity coupling approach and for the velocity-jump
  solver. */

  nv = 2;

  /**
  We change the dimension of the domain, the surface tension coefficient, and
  the coordinates of the center of the bubble. */

  f.sigma = 0.001;

  /**
  We must associate the height function field with the VOF fraction, so that it
  is used by the relevant functions (for curvature calculation). */

  f.height = h;

  /**
  We create the grid and start the simulation. */

  init_grid (1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {
  fraction (f, -circle(x,y,R0));

  /**
  We initialize the temperature field. The liquid phase temperature is set to
  the Scriven solution. */

  scalar TL = liq->T, TG = gas->T;
  foreach() {
    double r = sqrt (sq (x) + sq (y));
    TL[] = f[]*tempsol (r, R0);
    TG[] = (1. - f[])*TG0;
    T[]  = TL[] + TG[];
  }

  /**
  We set the boundary conditions for the liquid and gas phase temperature
  fields, which are those that are actually resolved by the phase change model.
  The one-field temperature `T` serves only for post-processing. */

  copy_bcs ({TL,TG}, T);
}

/**
We refine the domain according to the interface and the temperature field. */

#if TREE
event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet_leave_interface ({T,u}, {f},
      (double[]){1e-2,1e-2,uemax,uemax,uemax,1e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Logger

We output the total bubble volume in time (for testing). */

event logger (t += 0.1) {
  double bubblevol = 0.;
  foreach(reduction(+:bubblevol))
    bubblevol += (1. - f[])*dv();
  fprintf (stderr, "%d %.1f %.4f\n", i, t, bubblevol);
}

/**
### Movie

We write the animation with the evolution of the temperature field, the
gas-liquid interface and the grid refinement. */

event movie (t += 0.01, t <= 5) {
  clear();
  view (theta=0., phi=0., psi=-pi/2., width=1600.,
      tx = 0., ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = Tsat, max = Tbulk);
  mirror ({0.,1}) {
    draw_vof ("f", lw = 1.5);
    cells();
  }
  save ("movie.mp4");
}

