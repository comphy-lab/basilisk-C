/**
# Isothermal Evaporation of a Binary Droplet in Forced Convection

A binary liquid droplet, made of two components with the same
properties but with different volatilies evaporates in forced
convective conditions. The droplet is initially placed on the left
side of the domain. An inlet gas flowrate is imposed on the left
boundary, such that Re=160 with the initial diameter of the droplet.

The animation shows the evaporation of the liquid droplet, plotting
the mass fraction of the light component. The inlet velocity
transports the mass fraction toward the right boundary. The Reynolds
number selected for this simulation leads to the formation of
Von-Karman streets that can be visualized from the transport of the
chemical species mass fraction in gas phase.

![Evolution of the mass fraction fields of the light component](forcedbi/movie.mp4)
*/

/**
## Phase Change Setup

We let the default settings of the evaporation model: the
Stefan flow is shifted toward the liquid phase, and the
consistent phase is the liquid, which is advected using the
extended velocity. The multicomponent model requires the
number of gas and liquid species to be set as compiler
variables. We don't need to solve the temperature field
because the vapor pressure is set to a constant value,
different for each chemical species. */

#define NGS 3
#define NLS 2
#define FILTERED 1

/**
## Simulation Setup

We use the centered solver with the evaporation source term
in the projection step. The extended velocity is obtained
from the doubled pressure-velocity coupling. We use the
evaporation model together with the multiomponent phase
change mechanism. */

#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "multicomponent.h"
#include "balances.h"
#include "view.h"

/**
### Data for multicomponent model

We define the lists with the names of the chemical species
in gas and in liquid phase. The initial mass fractions are
defined for each component, as well as the diffusion
coefficients and the thermodynamic equilibrium constant.
Since we set the vapor pressure to a constant value, we
don't need to solve the temperature field. */

char* gas_species[NGS] = {"A", "B", "C"};
char* liq_species[NLS] = {"A", "B"};
char* inert_species[1] = {"C"};

double gas_start[NGS] = {0., 0., 1.0};
double liq_start[NLS] = {0.5, 0.5};

double inDmix1[NLS] = {1.4e-7, 1.4e-7};
double inDmix2[NGS] = {1.25e-5, 1.25e-5, 1.25e-5};

double inKeq[NLS] = {0.8, 0.4};

/**
### Boundary conditions

We set the inlet BCs on the left boundary, and outflow boundary
conditions on the right. Symmetry elsewhere. */

double vin = 1.424;
u.n[left] = dirichlet (vin);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);
uext.n[left] = dirichlet (vin);
uext.t[left] = dirichlet (0.);
pext[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation. */

int maxlevel, minlevel = 5;
double D0 = 0.4e-3, R0, effective_radius0;

int main (void) {
  /**
  We set the material properties of the fluids. */

  rho1 = 800., rho2 = 5.;
  mu1 = 1.138e-3, mu2 = 1.78e-5;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  R0 = 0.5*D0, L0 = 12.*(6.*R0);

  /**
  We change the surface tension coefficient. */

  f.sigma = 0.073;

  /**
  We run the simulation at different levels of
  refinement. */

  for (maxlevel = 9; maxlevel <= 9; maxlevel++) {
    CFL = 0.1;
    init_grid (1 << (maxlevel-3));
    run();
  }
}

#define circle(x, y, R) (sq(R) - sq(x - L0/6.) - sq(y - L0/2.))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

event init (i = 0) {
  refine (circle (x, y, 2.*R0) > 0. && level < maxlevel);
  fraction (f, circle (x, y, R0));
  effective_radius0 = sqrt (1./pi*statsf(f).sum);

#ifdef BALANCES
  mb.liq_species = liq_species;
  mb.gas_species = gas_species;
  mb.YLList = YLList;
  mb.YGList = YGList;
  mb.mEvapList = mEvapList;
  mb.liq_start = liq_start;
  mb.gas_start = gas_start;
  mb.rho1 = rho1;
  mb.rho2 = rho2;
  mb.inDmix1 = inDmix1;
  mb.inDmix2 = inDmix2;
  mb.maxlevel = maxlevel;
  mb.boundaries = true;
#endif
}

/**
We adapt the grid according to the mass fractions of the
species A and B, the velocity and the interface position. */

#if TREE
event adapt (i++) {
  scalar YA = YList[0], YB = YList[1];
  adapt_wavelet_leave_interface ({YA,YB,u.x,u.y}, {f},
      (double[]){1.e-3,1.e-3,1.e-3,1.e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Output Files

We write on a file the squared diameter decay and the dimensionless
time. */

event output_data (t += 5e-5) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double tad = t*inDmix2[0]/sq (2.*effective_radius0);
  double effective_radius = sqrt (statsf(f).sum/pi);
  double d_over_d0 = effective_radius/effective_radius0;
  double d_over_d02 = sq (d_over_d0);
  fprintf (fp, "%g %g %g %g\n", t, tad, d_over_d0, d_over_d02);
  fflush (fp);
}

/**
### Movie

We write the animation with the evolution of the light chemical
species mass fraction, and the interface position. */

event movie (t += 0.000125; t <= 0.03) {
  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f");
  squares ("A", min = 0., max = 0.5, linear = true);
  save ("movie.mp4");
}

/**
## Results

~~~gnuplot Squared Diameter Decay
reset
set xlabel "t [s/mm^2]"
set ylabel "(D/D_0)^2 [-]"

set size square
set key top right
set grid

plot "OutputData-9" u 2:4 w l lw 2 t "LEVEL 9"
~~~

~~~gnuplot Liquid Phase Mass Conservation
reset
set xlabel "t [s]"
set ylabel "(m_L - m_L^0) [kg]"
set key bottom left
set size square
set grid

plot "balances-9" every 500 u 1:10 w p ps 1.2 lc 1 title "Evaporated Mass Species A", \
     "balances-9" every 500 u 1:11 w p ps 1.2 lc 2 title "Evaporated Mass Species B", \
     "balances-9" every 500 u 1:4  w p ps 1.2 lc 3 title "Evaporated Mass Total", \
     "balances-9" u 1:(-$5) w l lw 2 lc 1 title "Variation Mass Species A", \
     "balances-9" u 1:(-$6) w l lw 2 lc 2 title "Variation Mass Species B", \
     "balances-9" u 1:(-$2) w l lw 2 lc 3 title "Variation Mass Total"
~~~

*/
