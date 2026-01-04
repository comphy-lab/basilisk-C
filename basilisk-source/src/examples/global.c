/**
# The Global Oceanic Circulation

We apply the generic ocean model to simulate the global oceanic
circulation. This is a "minimal model" similar to that used for the
[North Atlantic](gulf-stream.c).

The setup is close to that used for the global ocean model described in 
[Hurlburt & Hogan, 2000](#hurlburt2000) but uses the [layered
solver](/src/layered/hydro.h) described in [Popinet,
2020](/Bibliography#popinet2020).

See [Hurlburt & Hogan, 2000](#hurlburt2000) for details but the main
characteristics of the setup are:

* 6 isopycnal layers (see Table 2 of H&H, 2000 for the reference densities, thicknesses etc.)
* Wind stress in the top layer given by the monthly climatology of Hellerman & Rosenstein, 1983
* "Compressed bathymetry" as in H&H, 2000
* Quadratic bottom friction (Cb = 2 x 10^-3^), Laplacian horizontal viscosity (10 m^2^/s)

Note that this setup does not include the fluxes prescribed in H&H,
2000 on the northern boundary of the Atlantic. As a consequence the
AMOC and Gulf Stream are not very realistic in this global setup. It
should not be difficult to add these fluxes, but they are poorly
described in H&H, 2000.

![Animation of the norm of the surface velocity. Red is for values
 larger than 1.5 m/s. The spatial resolution is < 1/6
 degree.](global/2048/nu-short.mp4)(width=100%)

![Convergence with time of the average SSH](global/2048/etam.mp4)(width=100%)

# Setup

We use the generic ocean model with 6 isopycnal layers. The
documentation of the generic model needs to be consulted for
details. */

#define NL 6
#include "ocean.h"

/**
## Isopycnal layers

The setup uses six isopycnal layers with the initial thicknesses
specified by the `dh` array below. The corresponding relative density
differences are given by the `drho` array and correspond to those in
Table 2 of [Hurlburt & Hogan, 2000](#hurlburt2000). */

double * dh   = (double [NL]){ 0., 525, 375, 260, 185, 155 };
double * drho = (double [NL]){ 2.52/rho0, 2.28/rho0, 2.05/rho0, 1.78/rho0, 1.34/rho0, 0. };

/**
## Diapycnal entrainment

See [entrainment.h](/src/layered/entrainment.h) for explanations. The
mininum and maximum layer thicknesses are given by `hmin` and `hmax`
(in meters here). The average entrainement velocity is set to 0.07
cm/s (see Table 2 of H&H, 2000) and the coefficient of additional
interfacial friction associated with entrainment is zero. */

double * hmin = (double [NL]){ 40, 40, 40, 40, 40, 50 };
double * hmax = (double [NL]){ HUGE, HUGE, HUGE, HUGE, HUGE, HUGE };
double omr = 0.05e-2, Cm = 0.;

/**
## main() 

The simulation can be run in OpenMP/serial/GPU or with MPI. With MPI
the number of processes must be of the form 2^(i+1) i.e.: 8, 32, 128,
512, 2048 etc. (see [Tips](/src/Tips#non-cubic-domains) for
explanations). This can be done using the
[Makefile](/Tutorial#using-makefiles) with:

~~~bash
CC='mpicc -D_MPI=8' CFLAGS=-disable-dimensions make global.tst
~~~

it can also be run on GPUs using e.g.

~~~bash
OMP_NUM_THREADS=8 CFLAGS=-DSHOW make global.gpu.tst
~~~

Command-line parameters can be passed to the code to change the
spatial resolution and/or the timestep. */

int main (int argc, char * argv[])
{

  /**
  The domain is rectangular with a 2x1 aspect ratio. It is periodic in
  the longitude direction and spans [0:360]x[-90:90] degrees but is
  bounded to [-75:75] degrees in latitude using an artifical
  coastline. */

  dimensions (2, 1);
  size (360);
  origin (0, - 90);
  periodic (right);
  maxlat = 75;
  
  /**
  The Earth radius, here in meters, sets the length unit. */
  
  Radius = 6371220. [1];

  /**
  The acceleration of gravity sets the time unit (seconds here). */
  
  G = 9.8;

  /**
  The default resolution is 512 (longitude) x 256 (latitude)
  i.e. 180/256 $\approx$ 0.7 degree. */
  
  if (argc > 1)
    N = atoi(argv[1]);
  else
    N = 512;

  /**
  The default timestep is 600 seconds. */
  
  DT = 600;
  if (argc > 2)
    DT = atof(argv[2]);

  /**
  The default starting time for averaging (i.e. spinup time) is set to
  10 years (it can be overloaded with the third command line
  argument). */
  
  if (argc > 3)
    tspinup = atof(argv[3]);
  else
    tspinup = 10.*year;

  /**
  The number of layers is set to NL (six). */
  
  nl = NL;

  /**
  The "implicitness parameter" of the [implicit barotropic
  solver](/src/layered/implicit.h) is set to 0.55 rather than the
  default 0.5 ("Crank-Nicholson"). This causes numerical damping of
  the barotropic mode and is necessary to prevent "basin resonances"
  at low resolution (N = 512). At higher resolutions, values closer to
  0.5 (e.g. 0.51) seem to work fine. At all resolutions the
  sensitivity of the results to this parameter is low (resonances
  excepted). */
  
  theta_H = 0.55;

  run();
}

/**
## Run times

The simulation can/should be run on parallel machines using something like:

~~~bash
../qcc -source -D_MPI=1 global.c
scp _global.c navier.lmm.jussieu.fr:global/
mpicc -Wall -std=c99 -D_XOPEN_SOURCE=700 -O2 _global.c -o global -L$HOME/lib -lkdt -lm
~~~

On an [RTX
4090](https://www.techpowerup.com/gpu-specs/geforce-rtx-4090.c3889)
runtimes are approx. 70 years-per-day for a resolution of 2048 x 1024
(< 1/6 degree), i.e. less than 2 nanoseconds per grid point x step.

Spinup takes at least 10 years and statistics (for e.g. the standard
deviations of SSH) at least 10 years.

# References

~~~bib
@article{hurlburt2000,
  title={Impact of 1/8 to 1/64 resolution on {G}ulf {S}tream model--data 
         comparisons in basin-scale subtropical {A}tlantic {O}cean models},
  author={Hurlburt, Harley E and Hogan, Patrick J},
  journal={Dynamics of Atmospheres and Oceans},
  volume={32},
  number={3-4},
  pages={283--329},
  year={2000},
  publisher={Elsevier},
  DOI={10.1016/S0377-0265(00)00050-6},
  PDF={https://apps.dtic.mil/sti/tr/pdf/ADA531039.pdf}
}

@article{shriver1997,
  title={The contribution of the global thermohaline circulation to
  the Pacific to Indian Ocean throughflow via Indonesia},
  author={Shriver, Jay F and Hurlburt, Harley E},
  journal={Journal of Geophysical Research: Oceans},
  volume={102},
  number={C3},
  pages={5491--5511},
  year={1997},
  DOI={10.1029/1999GL011105},
  PDF={https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/1999GL011105}
}
~~~
*/
