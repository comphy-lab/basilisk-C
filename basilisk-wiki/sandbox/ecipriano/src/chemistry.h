/**
# Chemical Reactions

The chemical reactions are solved using the OpenSMOKE++ libraries
and an ODE system of equations. Despite the general formulation
of the chemistry, so far it is solved just in the gas phase,
focusing on the combustion kinetics. Different types of reactions
may require the implementation of new reactor systems 
([reactors.h](reactors.h)).
*/

#include "opensmoke.h"
#include "reactors.h"

/**
We initialize the OpenSMOKE++ ODE solver, which is the solver
that we are going to use for the chemistry integration. This
solver is specifically conceived for stiff reactive systems,
and it will automatically sub-step during the whole time step
of the simulation.
*/

event init (i = 0) {
  OpenSMOKE_InitODESolver ();
}

/**
We deallocate the OpenSMOKE++ ODE solver. */

event cleanup (t = end) {
  OpenSMOKE_CleanODESolver ();
}

/**
## Chemistry Event

This event assigns the batch reactor function, and it calls
the ODE system integrator advancing the chemical species
mass fraction and temperature fields according to the presence
of chemical reactions.
*/

event chemistry (i++) {

  /**
  Set reactor function and number of equations accordingly. */

#ifdef SOLVE_TEMPERATURE
  odefunction batch = &batch_nonisothermal_constantpressure;
  unsigned int NEQ = OpenSMOKE_NumberOfSpecies() + 1;
#else
  odefunction batch = &batch_isothermal_constantpressure;
  unsigned int NEQ = OpenSMOKE_NumberOfSpecies();
#endif

  foreach() {
    if (f[] < F_ERR) {

      /**
      Gather the initial conditions for the ODE solver. */

#ifdef SOLVE_TEMPERATURE
      double y0ode[NGS+1];
#else
      double y0ode[NGS];
#endif
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        y0ode[jj] = YG[];
      }
#ifdef SOLVE_TEMPERATURE
      y0ode[NGS] = TG[];
#endif

      /**
      Set additional data to be passed to the ODE system. */

      UserDataODE data;
      data.P = Pref;
      data.T = 1200.;
#ifdef SOLVE_TEMPERATURE
      data.rho = rho2;
      data.cp = cp2;
#endif
      double sources[NEQ];
      data.sources = sources;

      /**
      Solve the ODE system using the OpenSMOKE native ODE solver,
      specifically conceived for stiff reactive systems. */

      OpenSMOKE_ODESolver (batch, NEQ, dt, y0ode, &data);

      /**
      Recover the results of the ODE system. */

      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        YG[] = y0ode[jj];

#ifdef VARPROP
        scalar DYDt2jj = DYDt2[jj];
        DYDt2jj[] += sources[jj]*cm[];
#endif
      }
#ifdef SOLVE_TEMPERATURE
      TG[] = y0ode[NGS];
# ifdef VARPROP
      DTDt2[] += sources[NGS]*cm[];
# endif
#endif
    }
  }
}

