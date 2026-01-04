/**
# Bubble Rising in a Superheated Environment

Expansion of a rising bubble in a superheated environment.  A spherical bubble
is initialized in an AXI domain, in normal gravity conditions.  The bubble is at
saturation temperature while the liquid environment is superheated.  The phase
change pheomena, promoted by the temperature gradient between the surface of the
bubble and the liquid environment leads to the expansion of the bubble. At the
same time, the droplet rises the liquid column due to the presence of the
gravity.

The simulation setup used here was inspired
by [Bures et al., 2020](#burevs2021direct), who proposed this test
case, and by [Long et al., 2024](#long2024edge) who simulated the
[same case using EBIT](/sandbox/tianlong/test/bubblerising.c).

![Evolution of the interface and temperature field](rising/movie.mp4)(width="500" height="500")
*/

/**
## Simulation Setup

We use the centered Navier--Stokes equations solver with volumetric source in
the projection step. The phase change is directly included using the boiling
module, which sets the best (default) configuration for boiling problems. Many
features of the phase change (boiling) model can be modified directly in this
file without changing the source code, using the phase change model object
`pcm`. In this simulation we use a one-field velocity approach, and we shift the
volume expansion term toward the liquid phase. This allows using the
Navier--Stokes conserving approach, which guarantees conservation of the total
momentum. Therefore, we can combine the approaches by [Gennari et al. 2022](#gennari2022)
and [Boyd et al., 2023](#boyd2023). */

#include "axi.h"
#include "navier-stokes/low-mach.h"
#define FILTERED 1
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "boiling.h"
#include "view.h"

/**
### Boundary conditions

Outflow boundary conditions for velocity and pressure are imposed on the top,
left, and right walls. The temperature on these boundaries is imposed to the
bulk value. */

u.n[top] = dirichlet (0.);
u.t[top] = dirichlet (0.);
p[top] = neumann (0.);

u.n[left] = dirichlet (0.);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

double Tsat, Tbulk;
T[top] = dirichlet (Tbulk);
T[right] = neumann (0.);
T[left] = neumann (0.);

/**
### Problem Data

We declare the maximum and minimum levels of refinement,
the $\lambda$ parameter, the initial radius of the droplet,
the growth constant, and additional post-processing
variables. */

int maxlevel;
double R0 = 210e-6;
double XC, YC;
double betaGrowth;
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
  We set the material properties of the two fluids. In addition to the classic
  Basilisk setup for density and viscosity, we need to define thermal
  properties, such as the thermal conductivity $\lambda$, the heat capacity
  $cp$, and the enthalpy of vaporization $\Delta h_{ev}$. */

  rho1 = 757.0, rho2 = 1.435;
  mu1 = 4.29e-4, mu2 = 1.04e-5;
  lambda1 = 0.154, lambda2 = 0.02;
  cp1 = 3000., cp2 = 1830;
  dhev = 9.63e5;

  /**
  The initial bubble temperature and the interface temperature are set to the
  saturation value. */

  Tbulk = 354.55, Tsat = 351.45, TIntVal = Tsat;
  TL0 = Tbulk, TG0 = Tsat;

  /**
  In order to use the conserving extension of the Navier--Stokes solver we use
  a single pressure-velocity coupling (default), and we activate the `consistent`
  option of the phase change model. */

  nv = 1;
  pcm.consistent = true;

  /**
  We change the dimension of the domain and the origin of the bubble. */

  L0 = 20e-3 [*];
  XC = 1e-3, YC = 0.;

  /**
  We set the surface tension coefficient and the gravitational acceleration. */

  f.sigma = 0.018;
  G.x = -9.81;

  /**
  We run the simulation. */

  for (maxlevel = 10; maxlevel <= 10; maxlevel++) {
    init_grid (1 << 9);
    run();
  }
}

/**
A small circular bubble is initialized on the bottom side of the domain,
corresponding with the axis of symmetry. */

#define circle(x, y, R) (sq(R) - sq(x - XC) - sq(y - YC))

event init (i = 0) {

  /**
  In order to avoid initializing all domain at the maximum level of refinement,
  we refine only the region around the bubble. However, this procedure changes
  the values of fields defined by default in the `phase` model using TL0 and
  TG0.  Therefore, phase fields (TL and TG in this case) will have to be
  re-initialized. */

#if TREE
  refine (circle (x,y,4.*R0) > 0. && level < maxlevel);
#endif
  fraction (f, -circle(x,y,R0));

  double alpha = lambda1/rho1/cp1;
  betaGrowth = 0.5*R0/sqrt (alpha*0.0056);

  /**
  We initialize the temperature field. The liquid phase temperature is set to
  the analytic value. */

  scalar TL = liq->T, TG = gas->T;
  foreach() {
    double r = sqrt( sq(x - XC) + sq(y - YC) );
    TL[] = f[]*tempsol (r, R0);
    TG[] = (1. - f[])*TG0;
    T[]  = TL[] + TG[];
  }

  /**
  We set the boundary conditions for the liquid and gas phase temperature
  fields, which are those that are actually resolved by the phase change model.
  The one-field temperature `T` serves only for post-processing. */

  copy_bcs ({TL,TG}, T);

  /**
  Using a flux limiter we avoid spurious oscillations stemming from the
  discretization of the advection term. */

  phase_set_gradient (liq, minmod2);
  phase_set_gradient (gas, minmod2);
}

/**
We refine the domain according to the interface position, the temperature field,
and the liquid velocity. Increasing the maximum level of refinement, we can
encounter convergence issue, which are resolved by setting a higher minimum
level of refinement (5 in this case). */

#if TREE
event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet_leave_interface ({T,u.x,u.y}, {f},
      (double[]){1e-2,uemax,uemax}, maxlevel, 5, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Output Files

We calculate the effective bubble radius and the coordinate of its centroid. */

vector pos[];

event output (i++) {
  scalar fg[];
  foreach()
    fg[] = 1. - f[];

  foreach() {
    foreach_dimension()
      pos.x[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = mycs (point, f), p;
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      coord o = {x, y, z};
      foreach_dimension()
        pos.x[] = (o.x + p.x*Delta);
    }
  }

  scalar fx[], ff[];
  foreach() {
    ff[] = 1. - f[];
    fx[] = ff[]*x;
  }
  double xb = statsf(fx).sum / statsf(ff).sum;

  foreach()
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      pos.x[] -= xb;

  scalar ux[];
  foreach()
    foreach_dimension()
      ux[] = u.x[]*(1. - f[]);
  double ur = statsf(ux).sum;

  double Dx = statsf(pos.y).max;
  double Dy = statsf(pos.x).max - statsf(pos.x).min;
  double Ra = 0.25*(Dx + Dy);
  double Pe = rho1*cp1*Ra*ur/lambda1;
  Ra /= 2.*betaGrowth*sqrt (lambda1/rho1/cp1);

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);

  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g %g %g %g %g\n", t, effective_radius, Dx, Dy, Ra, Pe);
  fflush (fp);
}

/**
### Logger

We output the total bubble volume in time (for testing). */

event logger (t += 0.005) {
  double bubblevol = 0.;
  foreach(reduction(+:bubblevol))
    bubblevol += (1. - f[])*dv();
  fprintf (stderr, "%d %.3g %.3g\n", i, t, bubblevol);
}

/**
### Movie

We write the animation with the evolution of the temperature field and the
gas-liquid interface. */

event movie (t += 0.0002; t <= 0.08) {
  scalar TT[];
  foreach()
    TT[] = T[] - Tsat;

  clear();
  view (theta=0., phi=0., psi=-pi/2.,
        tx = 0., ty = 0.);

  travelling (0., 0.01, fov = 4, ty = -0.1);
  travelling (0.01, 0.8, tx = 0., ty = -5.5);

  draw_vof ("f", lw = 1.5);
  isoline ("TT", val = 2.8, filled = 1, fc={1.,1.,1.});
  squares ("T", min = Tsat, max = Tbulk, linear=true);
  mirror ({0.,1.}) {
    draw_vof ("f", lw = 1.5);
    isoline ("TT", val = 2.8, filled = 1, fc={1.,1.,1.});
    squares ("T", min = Tsat, max = Tbulk, linear=true);
    cells();
  }
  if (i > 1)
    save ("movie.mp4");
}

/**
## References

~~~bib
@article{burevs2021direct,
  title={Direct numerical simulation of evaporation and condensation with the geometric VOF method and a sharp-interface phase-change model},
  author={Bure{\v{s}}, Lubom{\'\i}r and Sato, Yohei},
  journal={International Journal of Heat and Mass Transfer},
  volume={173},
  pages={121233},
  year={2021},
  publisher={Elsevier}
}

@article{long2024edge,
  title={An Edge-based Interface Tracking (EBIT) Method for Multiphase Flows with Phase Change},
  author={Long, Tian and Pan, Jieyun and Zaleski, St{\'e}phane},
  journal={arXiv preprint arXiv:2402.13677},
  year={2024}
}
~~~
*/
