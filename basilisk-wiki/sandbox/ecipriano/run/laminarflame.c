/**
# Axisymmetric H2/N2–Air Coflow Flame

This test case reproduces a laminar, axisymmetric hydrogen–air diffusion flame
with detailed chemistry and transport, following the setup discussed by [Toro et al. 2004](#toro2005combined).

We use the low-mach solver in order to include weak compressibility effects due
to temperature and composition changes. The combustion module utilizes the phase
model for the solving the advection, diffusion, and reaction step for the
multicomponent gas mixture.
*/

#include "axi.h"
#include "navier-stokes/low-mach.h"
#include "navier-stokes/perfs.h"
#if USE_OPENSMOKE
# include "opensmoke/properties.h"
# include "opensmoke/chemistry.h"
#elif USE_CANTERA
# include "cantera/properties.h"
# include "cantera/chemistry.h"
#endif
#include "combustion.h"
#include "gravity.h"
#include "spark.h"
#include "view.h"

/**
## Simulation setup

We gather data for problem, including the kinetics folder with the detailed
hydrogen mechanism, the geometry of the problem, the inlet velocity, and the
inlet and initial species compositon (mass fractions). */

#define KINFOLDER "laminarflame.yaml"
#define FUEL_NOZZLE 9e-3
#define AIR_NOZZLE 95e-3
#define LENGTH 150e-3
#define U_AIR 0.5
#define H2_IN 0.066637
#define O2_IN 0.231442

/**
The inlet velocity is imposed using a parabolic profile, more realistic than a
flat profile for a laminar flow inside a nozzle. The central nozzle (y < R)
contains the fuel, while the outward nozzle corresponds to the nitrogen/oxigen
mixture. */

#define parabolic 2.*U_AIR*(1. - sq(y)/sq(R))
double R = 0.5*FUEL_NOZZLE, Re = 0.5*AIR_NOZZLE;

u.n[left] = dirichlet (y <= R ? parabolic : U_AIR);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

/**
Additional useful data are the maximum (and minimum) levels of refinemt, a
useful boolean for restarting long simulations, and a field which contains the
heat provided by the electric spark which ignites the mixture. */

int maxlevel, minlevel = 2;
bool restored = false;
scalar qspark[];

int main (int argc, char ** argv) {

  /**
  We set the kinetics folders, with the gas phase kinetics and the names of the
  chemical species in the scheme. Those are mostly used to facilitate pre- and
  post-processing operations. */

  kinetics (KINFOLDER, &NS);
  gas_species = new_species_names (NS);

  /**
  We set the initial thermodynamic state of the two-phase system (SI units). */

  Pref = 101325., T0 = 300.;

  /**
  We change the dimensions of the domain, introduce gravity, and we decide the
  maximum CFL number. */

  L0 = LENGTH;
  G.x = -9.81;
  CFL_MAX = 0.1;

  /**
  We initialize the grid and run the simulation. */

  for (maxlevel = 7; maxlevel <= 7; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }

  /**
  We de-allocate memory used for the vectors of species names. */

  free_species_names (NS, gas_species), gas_species = NULL;
}

event init (i = 0) {
  if (!restore (file = "restart")) {
    foreach()
      u.x[] = U_AIR;
  }
  else
    restored = true;

  /**
  We set the initial thermodynamic state of the mixture: T, P and mass
  fractions.  If the simulation is NOT restored, we force the initialization of
  the thermo state, otherwise, it would be effective only if the fields are
  constant. This is a guard to prevent problems when refining the grid in the
  `init` event (i.e. when using `refine()`). */

  double x[NS];
  x[index_species ("H2")] = 0.;
  x[index_species ("O2")] = O2_IN;
  x[index_species ("N2")] = 1. - O2_IN;

  //pcm.divergence = false;

  ThermoState tsg;
  tsg.T = T0, tsg.P = Pref, tsg.x = x;

  phase_set_thermo_state (gas, &tsg, force = !restored);

  /**
  The only property that we need to set, and that remains constant throughout
  the simulation, is the vector with the molecular weight of each species. */

  double MWGs[NS];
  molecular_weights (NS, MWGs);
  phase_set_properties (gas, MWs = MWGs);

  /**
  We set the spark properties: coordinate position, diameter, initial time,
  total duration, value, fields to which it is applied. Consider that the
  default policy is 3: heat source in the energy equation via `qspark`. Phase
  properties are also required. */

#if SPARK
  spark.T = qspark;
  spark.position = (coord){1.5e-3, 4.5e-3};
  spark.diameter = 1.5e-3;
  spark.time = 0.;
  spark.duration = 0.005;
  spark.temperature = 1e7;
  spark.phase = gas;
#endif

  /**
  The boundary conditions for species and temperature are set here because they
  do not exist in the global scope (i.e. where we set the velocity bcs), and
  because after the `init` event they would be skipped if the simulation is
  restored. */

  scalar fuel = gas->YList[index_species ("H2")];
  scalar oxi = gas->YList[index_species ("O2")];
  scalar inert = gas->YList[index_species ("N2")];
  scalar T = gas->T;

  fuel[left] = dirichlet (y <= R ? H2_IN : 0.);
  oxi[left] = dirichlet (y <= R ? 0. : O2_IN);
  inert[left] = dirichlet (y <= R ? 1. - H2_IN : 1. - O2_IN);
  T[left] = dirichlet (T0);
}

/**
### Mesh adaptation

We adapt the grid according to the fuel mass fraction, the temperature, and the
velocity field. */

#if TREE
event adapt (i++) {
  scalar fuel = gas->YList[index_species ("H2")];
  adapt_wavelet ({fuel,gas->T,u.x,u.y},
      (double[]){1e-2,1e0,1e-1,1e-1}, maxlevel, minlevel);
  unrefine (x >= (0.9*L0));
}
#endif

/**
## Tracing

Check the distribution of the computational time among the different events and
functions. */

#if TRACE > 1
event profiling (i += 20) {
  char name[80];
  sprintf (name, "profiling-%d", maxlevel);
  static FILE * fp = fopen (name, "w");
  trace_print (fp, 1);
}
#endif

/**
## Post-processing

The following events are for post-processing. */

/**
### Profiles

We extract the axial and radial profiles of temperature and mass fractions. */

event profiles (t = end) {
  char name[80];
  sprintf (name, "AxialProfiles-%d", maxlevel);
  FILE * faxial = fopen (name, "w");
  fprintf (faxial, "x(1) T(2) H2(3) H2O(4) O2(5) N2(6)\n");

  {
    coord p;
    coord box[2] = {{0, 0}, {140e-3,0}};
    coord n = {100,1};
    foreach_region (p, box, n) {
      fprintf (faxial, "%g %g %g %g %g %g\n", p.x,
          interpolate (gas->T, p.x, p.y),
          interpolate (gas->YList[index_species ("H2")], p.x, p.y),
          interpolate (gas->YList[index_species ("H2O")], p.x, p.y),
          interpolate (gas->YList[index_species ("O2")], p.x, p.y),
          interpolate (gas->YList[index_species ("N2")], p.x, p.y));
    }
  }

  char name3mm[80], name10mm[80], name20mm[80], name30mm[80];
  sprintf (name3mm,  "RadialProfiles3mm-%d",  maxlevel);
  sprintf (name10mm, "RadialProfiles10mm-%d", maxlevel);
  sprintf (name20mm, "RadialProfiles20mm-%d", maxlevel);
  sprintf (name30mm, "RadialProfiles30mm-%d", maxlevel);

  FILE * fp3mm  = fopen (name3mm, "w");
  FILE * fp10mm = fopen (name10mm, "w");
  FILE * fp20mm = fopen (name20mm, "w");
  FILE * fp30mm = fopen (name30mm, "w");

  {
    coord p;
    coord box[2] = {{3e-3,0}, {3e-3,15e-3}};
    coord n = {1,100};
    foreach_region (p, box, n)
      fprintf (fp3mm, "%g %g %g %g %g\n", p.y,
          interpolate (gas->T, p.x, p.y),
          interpolate (gas->XList[index_species ("H2")], p.x, p.y),
          interpolate (gas->XList[index_species ("O2")], p.x, p.y),
          interpolate (gas->XList[index_species ("H2O")], p.x, p.y));
  }

  {
    coord p;
    coord box[2] = {{10e-3,0}, {10e-3,15e-3}};
    coord n = {1,100};
    foreach_region (p, box, n)
      fprintf (fp10mm, "%g %g %g %g %g\n", p.y,
          interpolate (gas->T, p.x, p.y),
          interpolate (gas->XList[index_species ("H2")], p.x, p.y),
          interpolate (gas->XList[index_species ("O2")], p.x, p.y),
          interpolate (gas->XList[index_species ("H2O")], p.x, p.y));
  }

  {
    coord p;
    coord box[2] = {{20e-3,0}, {20e-3,15e-3}};
    coord n = {1,100};
    foreach_region (p, box, n)
      fprintf (fp20mm, "%g %g %g %g %g\n", p.y,
          interpolate (gas->T, p.x, p.y),
          interpolate (gas->XList[index_species ("H2")], p.x, p.y),
          interpolate (gas->XList[index_species ("O2")], p.x, p.y),
          interpolate (gas->XList[index_species ("H2O")], p.x, p.y));
  }

  {
    coord p;
    coord box[2] = {{30e-3,0}, {30e-3,15e-3}};
    coord n = {1,100};
    foreach_region (p, box, n)
      fprintf (fp30mm, "%g %g %g %g %g\n", p.y,
          interpolate (gas->T, p.x, p.y),
          interpolate (gas->XList[index_species ("H2")], p.x, p.y),
          interpolate (gas->XList[index_species ("O2")], p.x, p.y),
          interpolate (gas->XList[index_species ("H2O")], p.x, p.y));
  }

  fclose (faxial);
  fclose (fp3mm);
  fclose (fp10mm);
  fclose (fp20mm);
  fclose (fp30mm);

  scalar fuel = gas->YList[index_species ("H2")];
  scalar oxi = gas->YList[index_species ("O2")];
  scalar h2o = gas->YList[index_species ("H2O")];
  scalar T = gas->T;

  char name2[80];
  sprintf (name2, "Maps-%d", maxlevel);
  FILE * fpm = fopen (name2, "w");

  output_field ({fuel,oxi,h2o,T,u.x,u.y}, fp = fpm, linear = true,
      box = {{0,0}, {0.85*LENGTH,50e-3}});
}

/**
### Movie

Evolution of the temperature and the mass fractions of O2 and H2O in time. */

event movie (t += 0.01) {
  clear();
  view (tx = -0.5);
  squares ("T", min = 300., max = 1960.);
  mirror ({0,-1}) {
    cells();
  }
  save ("temperature.mp4");

  clear();
  view (tx = -0.5);
  squares ("YO2", min = 0., max = 0.232);
  mirror ({0,-1}) {
    squares ("YH2O", min = 0., max = 0.173);
  }
  save ("species.mp4");
}

event end (t = 0.3);

/**
## Results

We plot maps, radial, and axial profiles for species and temperatures. The comparison with experimental data shows displacements which can be due to:
(i) low level of refinement; (ii) no radiation if using Cantera; (iii) 
neglected Soret effect.

~~~gnuplot Temperature map
LEVEL = 7

set size ratio -1
unset key 
unset xtics
unset ytics
unset colorbox
set pm3d
set pm3d map interpolate 3,3
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
                          0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, \
                      0.625 1 0.9333 0, 0.75 1 0.4392 0, \
                      0.875 0.9333 0 0, 1 0.498 0 0 ) 

set xlabel "x"
set ylabel "y"
set title "temperature [K]"
set colorbox

splot "Maps-".LEVEL u 1:2:6
~~~

~~~gnuplot H2 map
set title "mass fraction H2 [-]"
splot "Maps-".LEVEL u 1:2:3
~~~

~~~gnuplot O2 map
set title "mass fraction O2 [-]"
splot "Maps-".LEVEL u 1:2:4
~~~

~~~gnuplot H2O map
set title "mass fraction H2O [-]"
splot "Maps-".LEVEL u 1:2:5
~~~

~~~gnuplot u.x map
set title "u.x [m/s]"
splot "Maps-".LEVEL u 1:2:7
~~~

~~~gnuplot Axial temperature profiles
reset
set grid
set xlabel "axial distance [mm]"
set ylabel "temperature [K]"

plot "AxialProfiles-".LEVEL u ($1*1e3):2 w l dt 1 lc -1 t "Temperature", \
     "../data/toro2004/50cms/axis-T.exp" u 1:2 w p lc -1 t "Toro et al., 2004 - CARS", \
     "../data/toro2004/50cms/axis-T.exp" u 3:4 w p lc -1 t "Toro et al., 2004 - Raman"
~~~

~~~gnuplot Radial temperature profiles at x = 3 mm
set grid
set xlabel "radial distance [mm]"
set ylabel "temperature [K]"
set xr[-15:15]
set yr[200:2200]

set y2tics
set y2r[0:1]
set y2label "mass fractions [-]"

plot "RadialProfiles3mm-".LEVEL u ($1*1e3):2 w l dt 1 lc -1 t "Temperature", \
     "RadialProfiles3mm-".LEVEL u (-$1*1e3):2 w l dt 1 lc -1 notitle, \
     "RadialProfiles3mm-".LEVEL u ($1*1e3):3 w l dt 1 lc 1 t "H2" axis x1y2, \
     "RadialProfiles3mm-".LEVEL u (-$1*1e3):3 w l dt 1 lc 1 notitle axis x1y2, \
     "RadialProfiles3mm-".LEVEL u ($1*1e3):4 w l dt 1 lc 2 t "O2" axis x1y2, \
     "RadialProfiles3mm-".LEVEL u (-$1*1e3):4 w l dt 1 lc 2 notitle axis x1y2, \
     "RadialProfiles3mm-".LEVEL u ($1*1e3):5 w l dt 1 lc 3 t "H2O" axis x1y2, \
     "RadialProfiles3mm-".LEVEL u (-$1*1e3):5 w l dt 1 lc 3 notitle axis x1y2, \
     "../data/toro2004/50cms/radial-3mm-T.exp" u 1:2 w p lc -1 t "Toro et al., 2004", \
     "../data/toro2004/50cms/radial-3mm-H2.exp" w p lc 1 axis x1y2 notitle, \
     "../data/toro2004/50cms/radial-3mm-O2.exp" w p lc 2 axis x1y2 notitle, \
     "../data/toro2004/50cms/radial-3mm-H2O.exp" w p lc 3 axis x1y2 notitle
~~~

~~~gnuplot Radial temperature profiles at x = 10 mm
set grid
set xlabel "radial distance [mm]"
set ylabel "temperature [K]"

set y2tics
set y2r[0:1]
set y2label "mass fractions [-]"

plot "RadialProfiles10mm-".LEVEL u ($1*1e3):2 w l dt 1 lc -1 t "Temperature", \
     "RadialProfiles10mm-".LEVEL u (-$1*1e3):2 w l dt 1 lc -1 notitle, \
     "RadialProfiles10mm-".LEVEL u ($1*1e3):3 w l dt 1 lc 1 t "H2" axis x1y2, \
     "RadialProfiles10mm-".LEVEL u (-$1*1e3):3 w l dt 1 lc 1 notitle axis x1y2, \
     "RadialProfiles10mm-".LEVEL u ($1*1e3):4 w l dt 1 lc 2 t "O2" axis x1y2, \
     "RadialProfiles10mm-".LEVEL u (-$1*1e3):4 w l dt 1 lc 2 notitle axis x1y2, \
     "RadialProfiles10mm-".LEVEL u ($1*1e3):5 w l dt 1 lc 3 t "H2O" axis x1y2, \
     "RadialProfiles10mm-".LEVEL u (-$1*1e3):5 w l dt 1 lc 3 notitle axis x1y2, \
     "../data/toro2004/50cms/radial-10mm-T.exp" u 1:2 w p lc -1 t "Toro et al., 2004 - CARS", \
     "../data/toro2004/50cms/radial-10mm-T.exp" u 3:4 w p lc -1 t "Toro et al., 2004 - Raman", \
     "../data/toro2004/50cms/radial-10mm-H2.exp" w p lc 1 axis x1y2 notitle, \
     "../data/toro2004/50cms/radial-10mm-O2.exp" w p lc 2 axis x1y2 notitle, \
     "../data/toro2004/50cms/radial-10mm-H2O.exp" w p lc 3 axis x1y2 notitle
~~~

~~~gnuplot Radial temperature profiles at x = 20 mm
set grid
set xlabel "radial distance [mm]"
set ylabel "temperature [K]"

set y2tics
set y2r[0:1]
set y2label "mass fractions [-]"

plot "RadialProfiles20mm-".LEVEL u ($1*1e3):2 w l dt 1 lc -1 t "Temperature", \
     "RadialProfiles20mm-".LEVEL u (-$1*1e3):2 w l dt 1 lc -1 notitle, \
     "RadialProfiles20mm-".LEVEL u ($1*1e3):3 w l dt 1 lc 1 t "H2" axis x1y2, \
     "RadialProfiles20mm-".LEVEL u (-$1*1e3):3 w l dt 1 lc 1 notitle axis x1y2, \
     "RadialProfiles20mm-".LEVEL u ($1*1e3):4 w l dt 1 lc 2 t "O2" axis x1y2, \
     "RadialProfiles20mm-".LEVEL u (-$1*1e3):4 w l dt 1 lc 2 notitle axis x1y2, \
     "RadialProfiles20mm-".LEVEL u ($1*1e3):5 w l dt 1 lc 3 t "H2O" axis x1y2, \
     "RadialProfiles20mm-".LEVEL u (-$1*1e3):5 w l dt 1 lc 3 notitle axis x1y2, \
     "../data/toro2004/50cms/radial-20mm-T.exp" u 1:2 w p lc -1 t "Toro et al., 2004 - CARS", \
     "../data/toro2004/50cms/radial-20mm-T.exp" u 3:4 w p lc -1 t "Toro et al., 2004 - Raman", \
     "../data/toro2004/50cms/radial-20mm-H2.exp" w p lc 1 axis x1y2 notitle, \
     "../data/toro2004/50cms/radial-20mm-O2.exp" w p lc 2 axis x1y2 notitle, \
     "../data/toro2004/50cms/radial-20mm-H2O.exp" w p lc 3 axis x1y2 notitle
~~~

~~~gnuplot Radial temperature profiles at x = 30 mm
set grid
set xlabel "radial distance [mm]"
set ylabel "temperature [K]"

set y2tics
set y2r[0:1]
set y2label "mass fractions [-]"

plot "RadialProfiles30mm-".LEVEL u ($1*1e3):2 w l dt 1 lc -1 t "Temperature", \
     "RadialProfiles30mm-".LEVEL u (-$1*1e3):2 w l dt 1 lc -1 notitle, \
     "RadialProfiles30mm-".LEVEL u ($1*1e3):3 w l dt 1 lc 1 t "H2" axis x1y2, \
     "RadialProfiles30mm-".LEVEL u (-$1*1e3):3 w l dt 1 lc 1 notitle axis x1y2, \
     "RadialProfiles30mm-".LEVEL u ($1*1e3):4 w l dt 1 lc 2 t "O2" axis x1y2, \
     "RadialProfiles30mm-".LEVEL u (-$1*1e3):4 w l dt 1 lc 2 notitle axis x1y2, \
     "RadialProfiles30mm-".LEVEL u ($1*1e3):5 w l dt 1 lc 3 t "H2O" axis x1y2, \
     "RadialProfiles30mm-".LEVEL u (-$1*1e3):5 w l dt 1 lc 3 notitle axis x1y2, \
     "../data/toro2004/50cms/radial-30mm-T.exp" u 1:2 w p lc -1 t "Toro et al., 2004 - CARS", \
     "../data/toro2004/50cms/radial-30mm-T.exp" u 3:4 w p lc -1 t "Toro et al., 2004 - Raman", \
     "../data/toro2004/50cms/radial-30mm-H2.exp" w p lc 1 axis x1y2 notitle, \
     "../data/toro2004/50cms/radial-30mm-O2.exp" w p lc 2 axis x1y2 notitle, \
     "../data/toro2004/50cms/radial-30mm-H2O.exp" w p lc 3 axis x1y2 notitle
~~~

## References

~~~bib
@article{toro2005combined,
  title={Combined experimental and computational study of laminar, axisymmetric hydrogen--air diffusion flames},
  author={Toro, VV and Mokhov, AV and Levinsky, HB and Smooke, MD},
  journal={Proceedings of the Combustion Institute},
  volume={30},
  number={1},
  pages={485--492},
  year={2005},
  publisher={Elsevier}
}
~~~
*/

