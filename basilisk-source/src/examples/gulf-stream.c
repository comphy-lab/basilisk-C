/**
# The Gulf Stream

We seek a "minimal setup" able to reproduce a realistic oceanic
circulation in the North Atlantic, with the Gulf Stream an obvious
dominant feature. The setup is as close as possible to that used by
[Hurlburt & Hogan, 2000](#hurlburt2000) but uses the [layered
solver](/src/layered/hydro.h) described in [Popinet,
2020](/Bibliography#popinet2020).

See [Hurlburt & Hogan, 2000](#hurlburt2000) for details but the main
characteristics of the setup are:

* 5 isopycnal layers (see Table 2 of H&H, 2000 for the reference
  densities, thicknesses etc.)
* Wind stress in the top layer given by the monthly climatology of
  Hellerman & Rosenstein, 1983
* "Compressed bathymetry" as in H&H, 2000
* Quadratic bottom friction (Cb = 2 x 10^-3^), Laplacian horizontal
  viscosity (10 m^2^/s)
* Atlantic Meridional Overturning Circulation (AMOC) driven by fluxes
  at the northern and southern boundaries (see Table 2 of H&H, 2000
  for the values of fluxes)

![Animation of the relative surface vorticity (approx. 2 years), min
 and max are $\pm$ 10^-4^ s^-1^. The spatial resolution is 1/24
 degree.](gulf-stream/2048/omega-short.mp4)(width=100%)

![Animation of the norm of the surface velocity. Red is for values
 larger than 1.5 m/s.](gulf-stream/2048/nu-short.mp4)(width=100%)

# Setup

We use the generic ocean model with 5 isopycnal layers. The
documentation of the generic model needs to be consulted for
details. */

#define NL 5
#include "ocean.h"

/**
## Isopycnal layers

The setup uses five isopycnal layers with an initial thickness of 250
meters for the top 4 layers (the `dh` array below). The corresponding
relative density differences are given by the `drho` array and
correspond to those in Table 2 of [Hurlburt & Hogan,
2000](#hurlburt2000). */

double * dh   = (double [NL]){ 0., 250, 250, 250, 250 };
double * drho = (double [NL]){ 2.13/rho0, 1.72/rho0, 1.41/rho0, 1.01/rho0, 0. };

/**
## Diapycnal entrainment

See [entrainment.h](/src/layered/entrainment.h) for explanations. The
mininum and maximum layer thicknesses are given by `hmin` and `hmax`
(in meters here). The average entrainement velocity is set to 0.07
cm/s (see Table 2 of H&H, 2000) and the coefficient of additional
interfacial friction associated with entrainment is zero. */

double * hmin = (double [NL]){ 40, 40, 40, 40, 50 };
double * hmax = (double [NL]){ HUGE, HUGE, HUGE, HUGE, HUGE };
double omr = 0.07e-2, Cm = 0.;

/**
## Boundary fluxes 

The [Atlantic Meridional Overturning Circulation
(AMOC)](https://en.wikipedia.org/wiki/Atlantic_meridional_overturning_circulation)
is driven primarily by heat (and salinity) fluxes at the surface of
the ocean, which are not included in this simplified model. Reasonably
realistic North Atlantic circulations can be obtained without them but
much more realistic results are obtained when including the simplified
representation proposed by [Hurlburt & Hogan,
2000](#hurlburt2000). See [bflux.h]() for details. */

#if 1
#include "bflux.h"
#endif

/**
## main() 

The simulation can be run in OpenMP/serial/GPU or with MPI. With MPI
the number of processes must be of the form 2^(i+1) i.e.: 8, 32, 128,
512, 2048 etc. (see [Tips](/src/Tips#non-cubic-domains) for
explanations). This can be done using the
[Makefile](/Tutorial#using-makefiles) with:

~~~bash
CC='mpicc -D_MPI=8' CFLAGS=-disable-dimensions make gulf-stream.tst
~~~

it can also be run on GPUs using e.g.

~~~bash
OMP_NUM_THREADS=8 CFLAGS=-DSHOW make gulf-stream.gpu.tst
~~~

Command-line parameters can be passed to the code to change the
spatial resolution and/or the timestep. */

int main (int argc, char * argv[])
{

  /**
  The domain is rectangular with a 2x1 aspect ratio. */  
  
  dimensions (2, 1);

  /**
  On a [spherical grid](/src/spherical.h), the sizes are given in
  degrees (of longitude). */
  
  size (84);

  /**
  The Earth radius, here in meters, sets the length unit. */
  
  Radius = 6371220 [1];

  /**
  The longitude and latitude of the lower-left corner. */
  
  origin (-98, 9);

  /**
  The acceleration of gravity sets the time unit (seconds here). */
  
  G = 9.8;

  /**
  The default resolution is 512 (longitude) x 256 (latitude)
  i.e. 42/256 $\approx$ 1/6 of a degree. */
  
  if (argc > 1)
    N = atoi(argv[1]);
  else
    N = 512;

  /**
  The default timestep is 600 seconds. Note that using a larger
  timestep (at low resolutions) can significantly affect the structure
  of the boundary current. */
  
  DT = 600 [0,1];
  if (argc > 2)
    DT = atof(argv[2]);

  /**
  We can change the spinup time using the third command-line parameter. */
  
  if (argc > 3)
    tspinup = atof(argv[3]);

  /**
  The number of layers is set to NL (five). */
  
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

The simulation can be run on parallel machines using something like:

~~~bash
../qcc -source -D_MPI=1 gulf-stream.c
scp _gulf-stream.c navier.lmm.jussieu.fr:gulf-stream/
mpicc -Wall -std=c99 -D_XOPEN_SOURCE=700 -O2 _gulf-stream.c -o gulf-stream -L$HOME/lib -lkdt -lm
~~~

On 128 cores of the "navier" cluster at d'Alembert, runtimes are of
the order of 89 simulated years per day (ypd) for a resolution of 512
x 256, 28 ypd for 1024 x 512 and 6 ypd for 2048 x 1024.

The simulation also runs fine [on GPUs](/src/grid/gpu/grid.h) using e.g.

~~~bash
OMP_NUM_THREADS=16 make gulf-stream.gpu.tst
~~~

On an [RTX 4090](https://www.techpowerup.com/gpu-specs/geforce-rtx-4090.c3889)
runtimes are approx. 48 ypd for a resolution of 2048 x 1024.

Spinup takes at least 5 years and statistics (for e.g. the standard
deviations of SSH) at least 10 years.

## Results

![Convergence with time of the average SSH](gulf-stream/2048/etam.mp4)(width=100%)

![Convergence with time of the standard deviation of SSH](gulf-stream/2048/etad.mp4)(width=100%)

All the reference generated files (log, perfs, *.mp4 etc.) are
accessible in the following folders (you will need to specify the
exact access path e.g. `../4096/eta.mp4`):

* [gulf-stream/512/]()
* [gulf-stream/1024/]()
* [gulf-stream/2048/]()
* [gulf-stream/4096/]()

## Todo

* Postprocessing (as in the [sandbox](/sandbox/popinet/gulf-stream.md)).
* Boundary fluxes cause large velocities which may unnecessarily
  restrict the timestep: smoother conditions may then improve runtimes.
* Real ("uncompressed") bathymetry runs but the Gulf Stream stays
  attached: this may be caused by a tuning of isopycnal layers, fluxes
  etc. which is specific to the compressed bathymetry.

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

@article{hurlburt2008,
  title={The {G}ulf {S}tream pathway and the impacts of the eddy-driven 
         abyssal circulation and the {D}eep {W}estern {B}oundary {C}urrent},
  author={Hurlburt, Harley E and Hogan, Patrick J},
  journal={Dynamics of Atmospheres and Oceans},
  volume={45},
  number={3-4},
  pages={71--101},
  year={2008},
  publisher={Elsevier},
  DOI={10.1016/j.dynatmoce.2008.06.002}
}
~~~
*/
