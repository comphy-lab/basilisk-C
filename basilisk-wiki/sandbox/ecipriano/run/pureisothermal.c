/**
# Isothermal Evaporation of a Pure Liquid Droplet

In this test case we want to simulate the evaporation of a
pure droplet in an isothermal environment. Even if the
droplet is pure, we use the multicomponent model to test
that, if the multiple species have identical properties,
they behave like a pure phase.

Since the evaporation is isothermal we don't need to use
*CLAPEYRON* or *ANTOINE*, but we just fix the value of the
thermodynamic equilibrium constant. For a pure droplet,
this is equivaluent to fixing the value of the chemical
species mass fraction on the gas-phase side of the
interface.

The animation shows the map of the mass fraction field in
gas phase, at different levels of refinement. The simulation
setup was borrowed from [Pathak et al., 2018](#pathak2018steady).

![Evolution of the mass fraction field](pureisothermal/movie.mp4)(height=500 width=500)
*/

/**
## Simulation Setup

We use the centered Navier--Stokes equations solver with volumetric source in
the projection step. The phase change is directly included using the evaporation
module, which sets the best (default) configuration for evaporation problems.
Many features of the phase change (evaporation) model can be modified directly
in this file without changing the source code, using the phase change model
object `pcm`. */

#include "navier-stokes/low-mach.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "view.h"

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right sides of the domain.
Boundary conditions for species must be set in the `init` event since those
fields are created in `defaults`. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
pf[top] = dirichlet_face (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
pf[right] = dirichlet_face (0.);

/**
### Simulation Data

We declare the maximum and minimum levels of refinement, the initial radius and
diameter, and the radius from the numerical simulation. */

int maxlevel, minlevel = 5;
double D0 = 0.4e-3, effective_radius0;

int main (void) {

  /**
  The number of gas and of liquid species are set in the `main()` function. */

  NGS = 5, NLS = 4;

  /**
  We set the material properties of the two fluids. In addition to the classic
  Basilisk setup for density and viscosity, we need to define species properties
  such as the diffusivity. The default thermodynamic pressure is the atmospheric
  value. To facilitate this setup we first set those properties which are common
  to all the species, and then we can refine the setup in the `init` event
  passing vectors to the `phase_set_properties()` function. */

  rho1 = 10.; rho2 = 1.;
  mu2 = 1.e-3; mu1 = 1.e-4;
  Dmix1 = 0., Dmix2 = 2.e-3;
  YIntVal = 0.667;

  /**
  We solve two different sets of Navier--Stokes equations according with the
  double pressure velocity coupling approach. The system is isothermal,
  therefore, the solution of the temperature equation is skipped. */

  nv = 2;
  pcm.isothermal = true;

  /**
  We change the dimension of the domain as a function of the initial diameter of
  the droplet. */

  L0 = 2.*D0;

  /**
  We change the surface tension coefficient. */

  f.sigma = 0.01;

  /**
  We run the simulation at different maximum levels of refinement. */

  for (maxlevel = 5; maxlevel <= 7; maxlevel++) {
    DT = 5e-8;
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the initial radius of the
droplet. We don't use the value D0 because for small errors of initialization
the squared diameter decay would not start from 1. */

event init (i = 0) {
  fraction (f, circle(x,y,0.5*D0));
  effective_radius0 = sqrt (4./pi*statsf(f).sum);

  /**
  The initial thermodynamics states of the two phases (i.e. the temperature,
  pressure, and composition) are defined and set. We force setting the thermo
  state only if the simulation was not restored. */

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = Pref, tsl.x = (double[]){0.25, 0.25, 0.25, 0.25};
  tsg.T = TG0, tsg.P = Pref, tsg.x = (double[]){0., 0., 0., 0., 1.};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  /**
  We set the boundary conditions for the mass fraction fields in liquid phase
  and for the inert. */

  for (scalar YG in gas->YList) {
    YG[top] = dirichlet (0.);
    YG[right] = dirichlet (0.);
  }
  scalar Inert = gas->YList[4];   // fixme: assuming species 4 is inert
  Inert[top] = dirichlet (1.);
  Inert[right] = dirichlet (1.);
}

/**
We adapt the grid according to the mass fraction of the evaporating species, and
the velocity field. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({Y,u}, {f},
      (double[]){1.e-3,1.e-2,1.e-2}, maxlevel, minlevel, 1);
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
  double d_over_d02 = sq (effective_radius / effective_radius0);
  double tad = t*Dmix2/sq(2.*effective_radius0);

  fprintf (fp, "%g %g %g\n", t, tad, d_over_d02);
  fflush (fp);
}

/**
### Logger

We output the total liquid volume in time (for testing). */

event logger (t += 1e-5) {
  fprintf (stderr, "%d %.1g %.3g\n", i, t, statsf (f).sum);
}

/**
### Mass Fraction Profiles

We write on a file the temperature and mass fraction profiles at different time
instants. */

event profiles (t = {1.03e-5, 6.03e-5, 1.40e-4}) {
  static FILE * fp = fopen ("Profiles", "w");

  static int pindex = 0;
  pindex = (pindex % 3) + 1;

  coord p;
  coord box[2] = {{0,0}, {L0,0}};
  coord n = {50,1};
  foreach_region (p, box, n)
    fprintf (fp, "maxlevel %d profile %d %g %g\n", maxlevel, pindex, p.x,
        interpolate (Y, p.x, p.y));
}

/**
### Movie

We write the animation with the evolution of the chemical species mass
fractions, the interface position and the temperature field. */

event movie (t += 0.5e-5, t <= 1.5e-4) {
  clear ();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("YG1", min = 0., max = 0.25*YIntVal, linear = true);
  save ("movie.mp4");
}

/**
## Results

The numerical results are compared with the results obtained
by [Pathak et al., 2018](#pathak2018steady) using a radially
symmetric model, integrated using an ODE solver.

~~~gnuplot Evolution of the squared diameter decay
reset
set xlabel "t [s]"
set ylabel "(D/D_0)^2"
set key top right
set grid

plot "../data/pathak-isothermal-unsteady-diam.csv" w p ps 2 t "Pathak et al., 2018", \
     "OutputData-5" u 1:3 w l lw 2 t "LEVEL 5", \
     "OutputData-6" u 1:3 w l lw 2 t "LEVEL 6", \
     "OutputData-7" u 1:3 w l lw 2 t "LEVEL 7"
~~~

~~~gnuplot Evolution of the chemical species mass fractions
reset
set xlabel "radius [m] x10^{6}"
set ylabel "Mass Fraction [-]"
set key top right
set grid

plot "../data/pathak-isothermal-unsteady-Yprofile-103e-5.csv" w p pt 8 lc 1 t "time = 1.03x10^{-5} s", \
     "../data/pathak-isothermal-unsteady-Yprofile-603e-5.csv" w p pt 8 lc 2 t "time = 6.03x10^{-5} s", \
     "../data/pathak-isothermal-unsteady-Yprofile-140e-4.csv" w p pt 8 lc 3 t "time = 1.40x10^{-5} s", \
     "<grep 'maxlevel 7 profile 1' Profiles" u ($5*1e+6):6 w l lw 2 lc 1 notitle, \
     "<grep 'maxlevel 7 profile 2' Profiles" u ($5*1e+6):6 w l lw 2 lc 2 notitle, \
     "<grep 'maxlevel 7 profile 3' Profiles" u ($5*1e+6):6 w l lw 2 lc 3 notitle
~~~

~~~bib
@article{pathak2018steady,
  title={Steady-state and transient solutions to drop evaporation in a finite domain: Alternative benchmarks to the d2 law},
  author={Pathak, Ashish and Raessi, Mehdi},
  journal={International Journal of Heat and Mass Transfer},
  volume={127},
  pages={1147--1158},
  year={2018},
  publisher={Elsevier}
}
~~~
*/
