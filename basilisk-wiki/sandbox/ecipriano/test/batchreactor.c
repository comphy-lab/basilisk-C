/**
# Batch Reactor - Non-Isothermal Constant Pressure

We test the solution of a non-isothermal, constant pressure,
batch reactor, using OpenSMOKE++ and the Basilisk interface.

The batch reactor is a perfectly stirred 0D system with
chemical reactions. The equations describing this system
are reported in [reactors.h](../src/reactors.h) and they
consists in an ODE system of equations solving chemical
species mass fractions and temperature.

We solve this system using the OpenSMOKE++ ODE solver,
specifically conceived for stiff reactive systems, and
we compare the results with a benchmark solution.

The system under investigation is the gas phase combustion
of hydrogen in air. The kinetic scheme features 32 species
and 173 reactions.
*/

/**
We define the number of gas species (I still can't get rid
of this variable, even if it's redundant since we have the
corresponding OpenSMOKE function). We aset the variable
properties solution in order to compute the density and the
heat capacity as a function of the thermodynamic state. */

#define NGS 32
#define VARPROP

/**
Most of the following lines are completely useless for the
solution of the problem, but we need them to avoid
compilation errors. */

attribute {
  scalar c;
}

char * gas_species[NGS] = {""};
scalar stefanflow[];

#include "fractions.h"
#include "common-evaporation.h"
#include "thermodynamics.h"
#include "opensmoke.h"
#include "reactors.h"

int main (void) {
  /**
  We initialize the OpenSMOKE pointers and we read the
  gas phase kinetic scheme. */

  OpenSMOKE_Init();
  OpenSMOKE_InitODESolver();
  OpenSMOKE_ReadKinetics ("/Users/ecipriano/OpenSMOKE/OpenSMOKEppTutorials/examples/OpenSMOKEpp_BatchReactor/01d-nonisothermal-constantpressure/kinetics-POLIMI_TOT_NOX_1412");

  /**
  We check that NGS is equal to the species found in the
  kinetic scheme (just a guard to avoid strange behavior).
  */

  assert (NGS == OpenSMOKE_NumberOfSpecies());

  /**
  We set the batch reactor equation function and the total
  number of equations/unknowns. */

  odefunction batch = &batch_nonisothermal_constantpressure;
  unsigned int NEQ = OpenSMOKE_NumberOfSpecies() + 1;

  /**
  We define an array with the initial conditions for every
  mass fraction value and for the temperature. */

  double y0unk[NEQ];
  for (int i=0; i<NEQ; i++)
    y0unk[i] = 0.;

  y0unk[OpenSMOKE_IndexOfSpecies ("N2")] = 7.452198e-01;
  y0unk[OpenSMOKE_IndexOfSpecies ("O2")] = 2.262686e-01;
  y0unk[OpenSMOKE_IndexOfSpecies ("H2")] = 2.851163e-02;
  y0unk[NGS] = 1000.;

  /**
  Data to be passed to the ODE system function are
  gathered in the following struct. We don't need to set
  the temperature because the reactor is non-isothermal;
  and we don't set rho and cp because we defined VARPROP.
  */

  UserDataODE data;
  data.P = 101325.;

  double sources[NEQ];
  data.sources = sources;

  /**
  We set the number of integration time steps, the total
  integration time and the time-steps (The ODE solver may
  internally perform more sub-steps if needed). */

  double nsteps = 10000.;
  double tend = 6.e-3;
  double deltat = tend/(double)nsteps;

  /**
  We prepare the output file header. */

  int counter = 3;
  fprintf (stderr, "time(1) T[K](2) ");
  for (int i=0; i<NGS; i++)
    fprintf (stderr, "%s(%d) ", OpenSMOKE_NamesOfSpecies (i), counter++);
  fprintf (stderr, "\n");

  /**
  We loop over time calling the ODE system integration. */

  for (double t=0; t<=tend; t+=deltat) {

    OpenSMOKE_ODESolver (batch, NEQ, deltat, y0unk, &data);

    /**
    Print species and temperature profiles. */

    fprintf (stderr, "%g %g ", t, y0unk[NGS]);
    for (int i=0; i<NGS; i++)
      fprintf (stderr, "%g ", y0unk[i]);
    fprintf (stderr, "\n");
  }

  /**
  Cleanup operations. */

  OpenSMOKE_Clean();
  OpenSMOKE_CleanODESolver();
}

/**
## Results

~~~gnuplot Temperature Profile
reset
set xr[0:0.0004]
set xlabel "t [s]"
set ylabel "Temperature [K]"
set key top left
set grid

plot "log" u 1:2 w l lw 2 lc 1 t "Results", \
     "batchreactor.data" u 1:5 w p lc 1 t "Reference"
~~~

~~~gnuplot Species Profiles
reset
set xr[0:0.0004]
set xlabel "t [s]"
set ylabel "Mass Fractions [-]"
set grid

plot "log" u 1:3 w l lw 2 lc 1 t "H2", \
     "log" u 1:4 w l lw 2 lc 2 t "O2", \
     "log" u 1:5 w l lw 2 lc 3 t "H2O", \
     "batchreactor.data" u 1:42 w p lc 1 notitle, \
     "batchreactor.data" u 1:43 w p lc 2 notitle, \
     "batchreactor.data" u 1:44 w p lc 3 notitle
~~~
*/
