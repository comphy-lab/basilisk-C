/**
# Isothermal Evaporation of a Static Binary Droplet

A binary liquid droplet is placed on the lower-left edge of the domain. The two
chemical species in liquid phase have the same physical properties, but
different volatility. The relative volatility between the heavy and the light
species is equal to 0.5. Therefore, we expect the light component to start to
evaporate first, increasing its mass fractions in gas phase and decreasing its
concentration in the liquid phase. The heavy component accumulates in the liquid
phase in response to the evaporation of the light species. The gas phase is
initially full of an inert compound, important in combustion simulations, which
always remains in gas phase.

![Evolution of the mass fraction fields of the light, heavy, and inert components, and the grid refinement](staticbi/movie.mp4)(height=600 width=600)
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

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

/**
### Simulation Data

We declare the maximum and minimum levels of refinement, the initial radius and
diameter, and the radius from the numerical simulation. */

int maxlevel, minlevel = 5;
double D0 = 0.4e-3, effective_radius0;

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

  rho1 = 10.; rho2 = 1.;
  mu1 = 1.e-4; mu2 = 1.e-5;
  Dmix1 = 4.e-6, Dmix2 = 8.e-5;

  /**
  We solve two different sets of Navier--Stokes equations according with the
  double pressure velocity coupling approach. The system is isothermal,
  therefore, the solution of the temperature equation is skipped. */

  nv = 2;
  pcm.isothermal = true;

  /**
  We change the dimension of the domain as a function of the initial diameter of
  the droplet. */

  D0 = 0.4e-3, L0 = 4.*D0;

  /**
  We change the surface tension coefficient. */

  f.sigma = 0.03;

  /**
  We run the simulation at three different levels of refinement. */

  for (maxlevel = 7; maxlevel <= 7; maxlevel++) {
    init_grid (1 << maxlevel);
    run ();
  }
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the initial radius of the
droplet. We don't use the value D0 because for small errors of initialization
the squared diameter decay would not start from 1. */

event init (i = 0) {
  fraction (f, circle (x, y, 0.5*D0));
  effective_radius0 = sqrt (4./pi*statsf(f).sum);

  /**
  The initial thermodynamics states of the two phases (i.e. the temperature,
  pressure, and composition) are defined and set. We force setting the thermo
  state only if the simulation was not restored. The composition is expressed in
  terms of mass fractions. The `YIntVals` vector sets the thermodynamic VLE
  equilibrium value for each chemical species. */

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = (double[]){0.5,0.5};
  tsg.T = TG0, tsg.P = Pref, tsg.x = (double[]){0.,0.,1};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  YIntVals[0] = 0.8, YIntVals[1] = 0.4;

  /**
  We overwrite the boundary conditions for the chemical species mass fractions
  in the *bcs* event. */

  for (scalar YG in gas->YList) {
    YG[top] = dirichlet (0.);
    YG[right] = dirichlet (0.);
  }
  scalar YI = gas->YList[NGS-1];
  YI[top] = dirichlet (1.);
  YI[right] = dirichlet (1.);
}

/**
We adapt the grid according to the mass fractions of the species A and B, the
velocity and the interface position. The vector Y is the sum of the evaporating
species mass fractions. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({Y,u.x,u.y}, {f},
      (double[]){1.e-4,1.e-3,1.e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Output Files

We write on a file the squared diameter decay and the dimensionless time. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double effective_radius = sqrt (4./pi*statsf(f).sum);
  double tad = t*Dmix2/sq (2.*effective_radius0);
  double d_over_d0  = effective_radius / effective_radius0;
  double d_over_d02 = sq (d_over_d0);

  fprintf (fp, "%g %g %g %g\n", t, tad, d_over_d0, d_over_d02);
  fflush (fp);
}

/**
### Logger

We output the total liquid volume in time (for testing). */

event logger (t += 1e-4) {
  fprintf (stderr, "%d %.3g %.3g\n", i, t, statsf (f).sum);
}

/**
### Movie

We write the animation with the evolution of the chemical species, the interface
position and the grid refinement. */

event movie (t += 2.e-5; t <= 0.005) {
  scalar YL1 = liq->YList[0], YL2 = liq->YList[1];
  scalar YG1 = gas->YList[0], YG2 = gas->YList[1];
  scalar YGI = gas->YList[2];

  scalar A[], B[], C[];
  foreach() {
    A[] = YL1[] + YG1[];
    B[] = YL2[] + YG2[];
    C[] = YGI[];
  }

  clear();
  draw_vof ("f", lw = 1.5);
  squares ("B", linear = true, min = 0., max = 0.56);
  mirror ({0,1}) {
    draw_vof ("f", lw = 1.5);
    squares ("C", linear = true, min = 0., max = 1.);
  }
  mirror ({1,0}) {
    draw_vof ("f", lw = 1.5);
    squares ("A", linear = true, min = 0., max = 0.5);
    mirror ({0,1}) {
      cells ();
      draw_vof ("f", lw = 1.5);
    }
  }
  save ("movie.mp4");
}

/**
## Results

~~~gnuplot Squared Diameter Decay
reset
set xlabel "t [s]"
set ylabel "(D/D_0)^2"
set key top right
set size square
set grid

plot "OutputData-7" u 2:4 w l lw 2 t "LEVEL 7"
~~~

The conservation tests compare the mass of the chemical species in
liquid phase with the total amount of the same species that
evaporates. If the global conservation is considered, the volume
fraction is used instead of the mass fraction field. See
[balances.h](../src/balances.h) for details.

~~~gnuplot Liquid Phase Mass Conservation
reset
set xlabel "t [s]"
set ylabel "(m_L - m_L^0) [kg]"
set key top right
set size square
set grid

plot "balances-7" every 500 u 1:10 w p ps 1.2 lc 1 title "Evaporated Mass Species A", \
     "balances-7" every 500 u 1:11 w p ps 1.2 lc 2 title "Evaporated Mass Species B", \
     "balances-7" every 500 u 1:4  w p ps 1.2 lc 3 title "Evaporated Mass Total", \
     "balances-7" u 1:(-$5) w l lw 2 lc 1 title "Variation Mass Species A", \
     "balances-7" u 1:(-$6) w l lw 2 lc 2 title "Variation Mass Species B", \
     "balances-7" u 1:(-$2) w l lw 2 lc 3 title "Variation Mass Total"
~~~

~~~gnuplot Gas Phase Mass Conservation
reset
set xlabel "t [s]"
set ylabel "(m_G - m_G^0) [kg]"
set key top left
set size square
set grid

plot "balances-7" every 500 u 1:(-$10) w p ps 1.2 lc 1 title "Evaporated Mass Species A", \
     "balances-7" every 500 u 1:(-$11) w p ps 1.2 lc 2 title "Evaporated Mass Species B", \
     "balances-7" every 500 u 1:(-$4)  w p ps 1.2 lc 3 title "Evaporated Mass Total", \
     "balances-7" u 1:7 w l lw 2 lc 1 title "Variation Mass Species A", \
     "balances-7" u 1:8 w l lw 2 lc 2 title "Variation Mass Species B", \
     "balances-7" u 1:3 w l lw 2 lc 3 title "Variation Mass Total"
~~~
*/

