/**
# Epstein Plesset

This test case was borrowed from
[sandbox/farsoiya](/sandbox/farsoiya/phase_change/epstein-plesset.c)
and it was adapted to my phase change model. It describes
the consumption of a liquid droplet in an isothermal environment,
where the evaporation is driven by a gradient of chemical species.

This test case considers mass transfer in pure diffusive conditions.
Therefore, the Stefan flow is not considered in the theretical
solution, obtained assuming quasi-static conditions. In this context,
the chemical species concentration profile is assumed to be equal to
the steady-state concentration profile, at any simulation time instants.
This approximation is justified when the diffusion is much faster than the interface regression velocity. Using this approximation, the theretical solution describing the evolution of the droplet radius in time was obtained by [Epstein and Plesset in 1950](#epstein1950stability):

$$
  \dfrac{dR}{dt} = -MW\dfrac{\mathcal{D} (\hat{c} - c_{bulk})}{\rho_g}
  \left(\dfrac{1}{R} + \dfrac{1}{\sqrt{\pi \mathcal{D} t}}\right)
$$

Even if [evaporation.h](../src/evaporation.h) is specifically
conceived for evaporation problems including the Stefan convection,
its formulation is general and it works also for diffusive
conditions. The initial setup of this simulation allows the expansion
term in the projection step to be null due to the density ratio equal
to 1, and the interface mass fraction is chosen to be sufficiently
small, so that the total vaporization rate tends to the pure
diffusive  conditions:

![Evolution of the Concentration Field](epsteinplesset/movie.mp4)
*/

/**
## Simulation Setup

We use the centered Navier--Stokes equations solver with volumetric source in
the projection step. The phase change is directly included using the evaporation
module, which sets the best (default) configuration for evaporation problems.
Many features of the phase change (evaporation) model can be modified directly
in this file without changing the source code, using the phase change model
object `pcm`. */

#include "axi.h"
#include "navier-stokes/low-mach.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "view.h"

/**
### Boundary Consitions

We let symmetry boundary conditions everywhere. The top and right mass fractions
are set to the bulk value. The absence of the Stefan flow in this simulation to
avoid open boundary conditions. */

Y[top] = dirichlet (0.);
Y[right] = dirichlet (0.);

/**
The original simulation setup is expressed in terms of chemical species
concentration, while the phase change model solves the mass fraction fields.
Therefore, we declare a concentration field used just for post-processing. */

scalar c[];

/**
### Problem Data

We declare the maximum and minimum levels of refinement, the initial droplet
radius and the radius from the numerical simulation. */

int maxlevel, minlevel = 5;
double R0 = 1., effective_radius0;

int main (void) {

  /**
  The number of gas and of liquid species are set in the `main()` function. */

  NGS = 2, NLS = 1;

  /**
  We set the material properties of the two fluids. In addition to the classic
  Basilisk setup for density and viscosity, we need to define species
  properties, such as the diffusivity $D$, and the molecular weights $MW$.
  */

  rho1 = 1., rho2 = 1.;
  mu2 = 1./20., mu1 = mu2/20.;
  Dmix1 = 0., Dmix2 = 1.;
  MW1 = 0.001, MW2 = MW1;

  /**
  We set the mass fraction of the evaporating species. The phase change model
  automatically adjust the value of the inert in order to close the mass
  fractions to 1. The interface is considered isothermal, and the thermodynamic
  VLE is a constant as well. */

  YG0 = 0., YL0 = 1.*MW1/rho1, YIntVal = YL0*0.8;

  /**
  We change the dimension of the domain, as a function of the initial droplet
  radius. */

  L0 = 10.*R0;

  /**
  Since there is no Stefan flow (the density ratio is 1), we can use a single
  velocity approach. The system is isothermal, therefore we skip the solution
  of the temperature equation. */

  nv = 1;
  pcm.isothermal = true;

  /**
  We change the surface tension coefficient of the droplet. The
  surface tension of the original test case is set to 0. A non-null
  surface tension is used here in order to reduce the time step,
  for the stability of the diffusion step, which includes an explicit
  source. */

  f.sigma = 0.1;

  /**
  We increase the grid size to speed up the simulation on the
  Basilisk wiki server. */

  for (maxlevel=8; maxlevel<=8; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

/**
In the initialization event, the volume fraction is initialized at the
lower-left corner of the domain, exploiting the spherical symmetry. We
initialize also he mass fractions in gas and in liquid phase, and we compute the
effective initial radius from the volume fraction field. */

event init (i = 0) {
  fraction (f, circle(x,y,R0));
  effective_radius0 = cbrt (3.*statsf(f).sum);

  /**
  We set the boundary conditions for the gas phase mass fraction, which is
  actually resolved by the phase change model. The one-field species `Y` serves
  only for post-processing. */

  scalar YG = gas->YList[0];
  copy_bcs ({YG}, Y);
}

/**
We refine the domain according to the interface and the concentration field. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({c,u.x,u.y}, {f},
      (double[]){1e-3,1e-3,1.e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Compute Concentration

After the solution of the tracer diffusion step, the concentration is
reconstructed from the mass fraction field. We do that in the *properties*
event, which is executed right after *tracer_diffusion*. */

event properties (i++) {
  scalar YG = gas->YList[0];
  foreach()
    c[] = YG[]*rho1/MW1 + f[]*1.;
}

/**
### Output File

We write on a file the droplet radius and the concentration of the chemical
species in a specific point of the domain, in time. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double effective_radius = cbrt (3.*statsf(f).sum);
  double radius_ratio = effective_radius / effective_radius0;

  fprintf (fp, "%g %g %g %g\n",
    t,  effective_radius, radius_ratio,
    interpolate(c, (1. + 0.2)*cos(M_PI/4), (1. + 0.2)*sin(M_PI/4) , 0) );
  fflush (fp);
}

/**
### Logger

We output the total liquid volume in time (for testing). */

event logger (t += 10) {
  fprintf (stderr, "%d %f %.3g\n", i, t, statsf (f).sum);
}

/**
### Movie

We write the animation with the evolution of the concentration field and the
gas-liquid interface. */

event movie (t += 10.; t <= 600.) {
  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("c", min = 0., max = 1., linear = true);
  save ("movie.mp4");
}

/**
## Results

~~~gnuplot Squared Diameter Decay
reset
set xlabel "t Dmix/D_0^2"
set ylabel "(D/D_0)^2"
set size square
set key top right
set grid

plot "../data/epsteinplesset.sol" u ($1/4):2 w p ps 1.4 title "Theretical", \
     "OutputData-8" u ($1/4):2 w l lw 2 title "LEVEL 8"
~~~

~~~gnuplot Evolution of the Concentration Field
reset
set xlabel "t Dmix/D_0^2"
set ylabel "Concentration [mol/m^3]"
set size square
set key top right
set grid

plot "../data/epsteinplesset.sol" u ($1/4):3 w p ps 1.4 title "Theoretical", \
     "OutputData-8" u ($1/4):4 w l lw 2 title "LEVEL 8"
~~~

## References

~~~bib
@article{epstein1950stability,
  title={On the stability of gas bubbles in liquid-gas solutions},
  author={Epstein, Paul S and Plesset, Milton S},
  journal={The Journal of Chemical Physics},
  volume={18},
  number={11},
  pages={1505--1509},
  year={1950},
  publisher={American Institute of Physics}
}
~~~
*/

