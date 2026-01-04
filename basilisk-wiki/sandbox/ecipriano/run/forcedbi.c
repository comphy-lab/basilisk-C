/**
# Isothermal Evaporation of a Binary Droplet in Forced Convection

A binary liquid droplet, made of two components with the same properties but
with different volatilies evaporates in forced convective conditions. The
droplet is initially placed on the left side of the domain. An inlet gas
flowrate is imposed on the left boundary, such that Re=160 with the initial
diameter of the droplet.

The animation shows the evaporation of the liquid droplet, plotting the mass
fraction of the light component. The inlet velocity transports the mass fraction
toward the right boundary. The Reynolds number selected for this simulation
leads to the formation of Von-Karman streets that can be visualized from the
transport of the chemical species mass fraction in gas phase.

![Evolution of the mass fraction fields of the light component](forcedbi/movie.mp4)
*/

/**
## Simulation Setup

We use the centered Navier--Stokes equations solver with volumetric source in
the projection step. The phase change is directly included using the evaporation
module, which sets the best (default) configuration for evaporation problems.
Many features of the phase change (evaporation) model can be modified directly
in this file without changing the source code, using the phase change model
object `pcm`. Exploiting the balances module we can automatically compute the
mass balances for each chemical species involved. */

#include "navier-stokes/low-mach.h"
#define FILTERED 1
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "balances/two-phase.h"
#include "view.h"

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right sides of the domain.
Boundary conditions for species must be set in the `init` event since those
fields are created in `defaults`. */

double vin = 1.424;
u.n[left] = dirichlet (vin);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

/**
### Simulation Data

We declare the maximum and minimum levels of refinement, the initial radius and
diameter, and the radius from the numerical simulation. */

int maxlevel, minlevel = 5;
double D0 = 0.4e-3, R0, effective_radius0;

int main (void) {

  /**
  The number of gas and of liquid species are set in the `main()` function. */

  NGS = 3, NLS = 2;

  /**
  We set the material properties of the two fluids. In addition to the classic
  Basilisk setup for density and viscosity, we need to define species properties
  such as the diffusivity. The default thermodynamic pressure is the atmospheric
  value. To facilitate this setup we first set those properties which are common
  to all the species, and then we can refine the setup in the `init` event
  passing vectors to the `phase_set_properties()` function. */

  rho1 = 800., rho2 = 5.;
  mu1 = 1.138e-3, mu2 = 1.78e-5;
  Dmix1 = 1.4e-7, Dmix2 = 1.25e-5;

  /**
  We solve two different sets of Navier--Stokes equations according with the
  double pressure velocity coupling approach. The system is isothermal,
  therefore, the solution of the temperature equation is skipped. */

  nv = 2;
  pcm.isothermal = true;

  /**
  We change the dimension of the domain as a function of the initial diameter of
  the droplet. */

  R0 = 0.5*D0, L0 = 12.*(6.*R0);

  /**
  We change the surface tension coefficient. */

  f.sigma = 0.073;

  /**
  We run the simulation at different levels of refinement. */

  for (maxlevel = 9; maxlevel <= 9; maxlevel++) {
    init_grid (1 << (maxlevel-3));
    run();
  }
}

#define circle(x, y, R) (sq(R) - sq(x - L0/6.) - sq(y - L0/2.))

/**
We initialize the volume fraction field and we compute the initial radius of the
droplet. We don't use the value D0 because for small errors of initialization
the squared diameter decay would not start from 1. */

event init (i = 0) {
  refine (circle (x, y, 2.*R0) > 0. && level < maxlevel);
  fraction (f, circle (x, y, R0));
  effective_radius0 = sqrt (1./pi*statsf(f).sum);

  /**
  The initial thermodynamics states of the two phases (i.e. the temperature,
  pressure, and composition) are defined and set. We force setting the thermo
  state only if the simulation was not restored. The composition is expressed in
  terms of mass fractions. The `YIntVals` vector sets the thermodynamic VLE
  equilibrium value for each chemical species. */

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = (double[]){0.5,0.5};
  tsg.T = TG0, tsg.P = Pref, tsg.x = (double[]){0.,0.,1};

  phase_set_thermo_state (liq, &tsl, force = true);
  phase_set_thermo_state (gas, &tsg, force = true);

  YIntVals[0] = 0.8, YIntVals[1] = 0.4;
}

/**
We adapt the grid according to the mass fractions of the species A and B, the
velocity and the interface position. The vector Y is the sum of the evaporating
species mass fractions. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({Y,u.x,u.y}, {f},
      (double[]){1.e-3,1.e-3,1.e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Output Files

We write on a file the squared diameter decay and the dimensionless time. */

event output_data (t += 5e-5) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double tad = t*Dmix2/sq (2.*effective_radius0);
  double effective_radius = sqrt (statsf(f).sum/pi);
  double d_over_d0 = effective_radius/effective_radius0;
  double d_over_d02 = sq (d_over_d0);
  fprintf (fp, "%g %g %g %g\n", t, tad, d_over_d0, d_over_d02);
  fflush (fp);
}

/**
### Logger

We output the total liquid volume in time (for testing). */

event logger (t += 1e-3) {
  fprintf (stderr, "%d %.4g %.4g\n", i, t, statsf (f).sum);
}

/**
### Movie

We write the animation with the evolution of the light chemical species mass
fraction, and the interface position. */

event movie (t += 0.000125; t <= 0.03) {
  scalar A[], AL = liq->YList[0], AG = gas->YList[0];
  foreach()
    A[] = AL[] + AG[];

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
