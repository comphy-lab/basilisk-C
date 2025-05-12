/**
# Rising of a Growing Bubble

Expansion of a rising bubble in a superheated environment.
A spherical bubble is initialized in an AXI domain,
in normal gravity conditions.  The bubble is at saturation
temperature while the liquid environment is superheated.
The phase change pheomena, promoted by the temperature
gradient between the surface of the bubble and the liquid
environment leads to the expansion of the bubble. At the
same time, the droplet rises the liquid column due to the
presence of the gravity.

The simulation setup used here was inspired
by [Tanguy et al., 2014](#tanguy2014benchmarks), and
is characterized by $\text{Ja}=3$. There is no analytic solution.

The animation shows the map of the temperature field, and
the evolution of the gas-liquid interface. Two different
situations are considered:

1. High surface tension coefficient: $\sigma=0.07 \text{ Nm}^{-1}$,
   which forces the bubble to remain spherical during the
   entire simulation ($\text{We} < 0.1$).

2. Low surface tension coefficient $\sigma=0.001 \text{ Nm}^{-1}$,
   which leads to the bubble deformation.

![Evolution of the interface and temperature field $\sigma=0.001$](bubblerisengrowth/movie1.mp4)(width="500" height="500")
![Evolution of the interface and temperature field $\sigma=0.07$](bubblerisengrowth/movie2.mp4)(width="500" height="500")
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
#define SOLVE_LIQONLY

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
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
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
imposed on the top, left, and right walls. The temperature
on these boundaries is imposed to the bulk value. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
uext.n[top] = neumann (0.);
uext.t[top] = neumann (0.);
pext[top] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);
uext.n[left] = neumann (0.);
uext.t[left] = neumann (0.);
pext[left] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);

TL[top] = dirichlet (Tbulk);
TL[right] = neumann (0.);
TL[left] = neumann (0.);

/**
### Problem Data

We declare the maximum and minimum levels of refinement,
the $\lambda$ parameter, the initial radius of the droplet,
the growth constant, and additional post-processing
variables. */

int maxlevel = 8, minlevel = 5, sim;
double R0 = 0.1e-3;
double XC, YC;
double betaGrowth = 3.32615013;
double effective_radius;

/**
### Initial Temperature

We set the initial temperature profile to the Scriven
solution. */

#include <gsl/gsl_integration.h>
#pragma autolink -lgsl -lgslcblas

double intfun (double x, void * params) {
  double beta = *(double *) params;
  return exp(-sq(beta)*(pow(1. - x, -2.) - 2.*(1. - rho2/rho1)*x - 1 ));
}

double tempsol (double r, double R) {
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  double result, error;
  double beta = betaGrowth;
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
  We set the material properties of the fluids. */

  rho1 = 958.; rho2 = 0.59;
  mu1 = 2.82e-4; mu2 = 1.23e-6;
  lambda1 = 0.6, lambda2 = 0.026;
  cp1 = 4216., cp2 = 2034.;
  dhev = 2.257e+6;

  /**
  The initial bubble temperature and the interface
  temperature are set to the saturation value. */

  Tbulk = 373.989, Tsat = 373., TIntVal = 373.;
  TL0 = Tbulk, TG0 = TIntVal;

  /**
  We change the dimension of the domain,
  the surface tension coefficient, and the coordinates
  of the center of the bubble. */

  L0 = 2.4e-3;
  XC = 0.15*L0, YC = 0.;

  /**
  We reduce the tolerance of the Poisson equation solver,
  and the maximum allowed time step. */

  TOLERANCE = 1.e-6;
  DT = 1.e-5;

  /**
  We add the gravity contribution using the
  reduced approach, which applied the gravity
  force just at the gas-liquid interface. */

  G.x = -9.81;

  /**
  We define a list with the surface tension coefficients
  used in the two different simulations. */

  double sigmas[2] = {0.001, 0.07};

  /**
  We set the surface tension and run the simulation. */

  for (sim=0; sim<=1; sim++) {
    f.sigma = sigmas[sim];
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x, y, R) (sq(R) - sq(x - XC) - sq(y - YC))

event init (i = 0) {
  fraction (f, -circle(x,y,R0));

  /**
  We initialize the temperature field. The liquid
  phase temperature is set to the analytic value. */

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
  double uemax = 1e-4;
  adapt_wavelet_leave_interface ({T,u.x,u.y}, {f},
      (double[]){1e-4,uemax,uemax}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Output Files

We reconstruct the effective bubble radius and we write it
on a file. */

event output (i++) {
  scalar fg[];
  foreach()
    fg[] = 1. - f[];

  effective_radius = pow(3./2.*statsf(fg).sum, 1./3.);

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);

  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t, effective_radius);
  fflush (fp);
}

/**
### Movie

We write the animation with the evolution of the
temperature field and the gas-liquid interface. */

event movie (t += 0.0002; t <= 0.02) {
  clear();
  view (theta=0., phi=0., psi=-pi/2.,
        tx = 0., ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = Tsat, max = Tbulk, linear=true);
  mirror ({0.,1.}) {
    draw_vof ("f", lw = 1.5);
    squares ("T", min = Tsat, max = Tbulk, linear=true);
  }
  if (sim == 0)
    save ("movie1.mp4");
  else
    save ("movie2.mp4");
}

/**
## References

~~~bib
@article{tanguy2014benchmarks,
  title={Benchmarks and numerical methods for the simulation of boiling flows},
  author={Tanguy, S{\'e}bastien and Sagan, Micha{\"e}l and Lalanne, Benjamin and Couderc, Fr{\'e}d{\'e}ric and Colin, Catherine},
  journal={Journal of Computational Physics},
  volume={264},
  pages={1--22},
  year={2014},
  publisher={Elsevier}
}
~~~
*/
