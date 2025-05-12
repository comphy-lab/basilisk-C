/**
# Expansion of a Nucleated Bubble

A spherical bubble is initialized on the lower-left
corner of and AXI domanin. Similarly to the Scriven
problem, the bubble is expanded due to the temperature
gradient between the interface of the bubble, which
remains at saturation temperature, and the superheated
liquid environment. However, in this case we consider
the left side of the boundary as a solid wall and we
impose a contact angle, which influences the bubble
growth process.

This is the simplest nucleate boiling configuration,
where no microlayer modelling and no heat exchange
with the solid wall is condered.

![Evolution of the interface, the temperature field, and the grid refinement](bubblecontact/movie.mp4)(height=400 width=900)
*/

/**
## Phase Change Setup

We move the interface using the velocity *uf*, with the
expansion term shifted toward the gas-phase. In this way
*uf* is divergence-free at the interface. The double
pressure velocity couping is used to obtain an extended
velocity, used to transport the gas phase tracers. */

#define INT_USE_UF
#define CONSISTENTPHASE2
#define SHIFT_TO_GAS
#define INIT_TEMP

/**
## Simulation Setup

We use the centered solver with the divergence source term,
and the extended velocity is obtained using the centered-doubled
procedure. The evaporation model is used in combination
with the temperature-gradient mechanism, which manages
the solution of the temperature field. */

#include "axi.h"
#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "temperature-gradient.h"
#include "view.h"

/**
We declare the variables required by the
temperature-gradient model. */

double lambda1, lambda2, cp1, cp2, dhev;
double TL0, TG0, TIntVal, Tsat, Tbulk;

/**
### Boundary conditions

Outflow boundary conditions for velocity and pressure are
imposed on the top and right walls. On the left wall,
no-slip conditions are imposed.*/

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
uext.n[top] = neumann (0.);
uext.t[top] = neumann (0.);
pext[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);

u.n[left] = dirichlet (0.);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);
uext.n[left] = dirichlet (0.);
uext.t[left] = dirichlet (0.);
pext[left] = neumann (0.);

TL[top] = dirichlet (Tbulk);
TL[right] = dirichlet (Tbulk);

/**
To set the contact angle, we allocate a [height-function
field](/src/heights.h) and set the contact angle boundary condition on
its tangential component. */

vector h[];
double theta0 = 40.;
h.t[left] = contact_angle (theta0*pi/180.);

/**
### Problem Data

We declare the maximum and minimum levels of refinement,
the $\lambda$ parameter, the initial radius of the droplet,
the growth constant, and additional post-processing
variables. */

int maxlevel = 7, minlevel = 3;
double R0 = 0.1;
double XC, YC;
double V0, tshift;
double betaGrowth = 0.7659010001953077;

/**
### Initial Temperature

We use gsl to compute the integral required to
obtain the Scriven temperature profile.
*/

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
  return Tbulk - 2.*sq(beta)*(rho2*(dhev + (cp1 - cp2)*(Tbulk - Tsat))/rho1/cp1)*result;
}

int main (void) {
  /**
  We set the material properties of the fluids (Ja=0.5). */

  rho1 = 2.5; rho2 = 0.25;
  mu1 = 7.e-3; mu2 = 7.e-4;
  lambda1 = 0.07, lambda2 = 0.007;
  cp1 = 2.5, cp2 = 1.;
  dhev = 100.;

  /**
  The initial bubble temperature and the interface
  temperature are set to the saturation value. */

  Tbulk = 3., Tsat = 1., TIntVal = 1.;
  TL0 = Tbulk, TG0 = TIntVal;

  /**
  We change the dimension of the domain,
  the surface tension coefficient, and the coordinates
  of the center of the bubble. */

  L0 = 1.;
  XC = 0., YC = 0.;
  f.sigma = 0.001;

  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */

  f.height = h;

  /**
  We create the grid and start the simulation. */

  init_grid (1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x - XC) - sq(y - YC))

event init (i = 0) {
  fraction (f, circle(x,y,R0));
  foreach()
    f[] = 1. - f[];

  /**
  We initialize the temperature field. The liquid
  phase temperature is set to the Scriven solution. */

  foreach() {
    double r = sqrt( sq(x - XC) + sq(y - YC) );
#ifdef INIT_TEMP
      TL[] = f[]*tempsol (r, R0);
#else
    TL[] = f[]*TL0;
#endif
    TG[] = (1. - f[])*TG0;
    T[]  = TL[] + TG[];
  }
}

/**
We refine the domain according to the interface and the
temperature field. */

#if TREE
event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet_leave_interface ({T,u}, {f},
      (double[]){1e-2,uemax,uemax,uemax,1e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Movie

We write the animation with the evolution of the
temperature field, the gas-liquid interface and
the grid refinement. */

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

