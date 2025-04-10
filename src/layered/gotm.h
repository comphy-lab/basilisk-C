/**
# General Ocean Turbulence Model interface for multilayer

This module interfaces the [multilayer solver](hydro.h) with
[GOTM](http://gotm.net).

Vertical diffusion of momentum, temperature and salinity is handled
by GOTM.

The temperature and salinity fields are undefined by default. If the
calling solver defines them, the corresponding GOTM vertical diffusion
routines will be called. */

scalar T = {-1}, S = {-1};

/**
The surface fluxes of momentum, heat and short-wave radiation are
provided as surface fields by the calling solver. They are zero by
default. */

// surface wind stress (eg. in Pa or N/m^2)
(const) vector airsea_tau[] = {0.};
// surface heat flux (eg. in W/m^2)
(const) scalar airsea_heat_flux[] = 0.;
// surface short-wave radiation flux (eg. in W/m^2)
(const) scalar airsea_swr_flux[] = 0.;

/**
## Implementation

It uses the relevant headers from the [C
interface](/src/gotm/common.h) to GOTM. */

#include "gotm/common.h"
#include <gotm/gotm/gotm.h>
// #include <gotm/gotm/diagnostics.h>
// #include <gotm/util/time.h>

#include <gotm/meanflow/meanflow.h>
// #include <gotm/meanflow/updategrid.h>
// #include <gotm/meanflow/wequation.h>
#include <gotm/meanflow/coriolis.h>
#include <gotm/meanflow/uequation.h>
#include <gotm/meanflow/vequation.h>
// #include <gotm/meanflow/extpressure.h>
// #include <gotm/meanflow/intpressure.h>
#include <gotm/meanflow/friction.h>
#include <gotm/meanflow/shear.h>
#include <gotm/meanflow/stratification.h>
#include <gotm/meanflow/salinity.h>
#include <gotm/meanflow/temperature.h>
#include <gotm/meanflow/convectiveadjustment.h>

#include <gotm/util/convert_fluxes.h>
#include <gotm/util/tridiagonal.h>
#include <gotm/util/eqstate.h>

#include <gotm/turbulence/turbulence.h>
#include <gotm/turbulence/kpp.h>

// #include <gotm/observations/observations.h>
#include <gotm/airsea/airsea.h>

/**
The corresponding GOTM libraries need to be linked. Note that the
`GOTM` environment variable need to point to the proper installation
directory. */

#pragma autolink							\
  -L$GOTM/build/ -lgotm	-lairsea -lmeanflow -lobservations -lturbulence \
  -linput -lutil -loutput_manager  `nf-config --flibs`			\
  -lgfortran

/**
The timestepping routine is based on its [Fortran
counterpart](https://github.com/gotm-model/code/blob/v5.2/src/gotm/gotm.F90#L372). */

static void gotm_step (long n)
{
  // prepare time and output
  //  time_update_time (&n);

  //     all observations/data

#if 0  
  input_do_input (&time_julianday, &time_secondsofday,
		  &nl, &meanflow_z);
  observations_get_all_obs (&time_julianday, &time_secondsofday,
			    &nl, &meanflow_z);
#endif
  
  //     external forcing
#if 0  
  if (airsea_calc_fluxes) {
    // fixme: check indices
    airsea_set_sst (&meanflow_t.a[nl-1]); // check this
    airsea_set_ssuv (&meanflow_u.a[nl-1],
		     &meanflow_v.a[nl-1]);
  }
#endif
  //  airsea_do_airsea (&time_julianday, &time_secondsofday);

  //     reset some quantities
  airsea_tx /= meanflow_rho_0;
  airsea_ty /= meanflow_rho_0;

  // sum fluxes as a diagnostic
  //  airsea_integrated_fluxes (&dt);

  //     meanflow integration starts
#if 0  
  realtype zeta = observations_get_zeta();
  meanflow_updategrid (&nl, &dt, &zeta);
  meanflow_wequation (&nl, &dt);
#endif

#if 0
  double omega = 2.*pi/86164.;
  double latitude = 50.; // 50.;
  meanflow_cori = 2.*omega*sin(2.*pi*latitude/360.);
  meanflow_coriolis (&nl, &dt);
#endif
  
  //     update velocity
  integer ext_press_mode = 99; // no external pressure gradient
  meanflow_uequation (&nl, &dt, &gotm_cnpar, &airsea_tx,
		      turbulence_num.a, turbulence_gamu.a,
		      &ext_press_mode);
  meanflow_vequation (&nl, &dt, &gotm_cnpar, &airsea_ty,
		      turbulence_num.a, turbulence_gamv.a,
		      &ext_press_mode);
#if 0
  meanflow_extpressure (&observations_ext_press_mode, &nl);
  meanflow_intpressure (&nl);
#endif
  meanflow_friction (&turbulence_kappa, &meanflow_avmolu,
		     &airsea_tx, &airsea_ty);

  //     update temperature and salinity
  if (S.i >= 0)
    meanflow_salinity (&nl, &dt, &gotm_cnpar,
		       turbulence_nus.a, turbulence_gams.a);
  
  if (T.i >= 0)
    meanflow_temperature (&nl, &dt, &gotm_cnpar, &airsea_i_0, &airsea_heat,
			  turbulence_nuh.a, turbulence_gamh.a, meanflow_rad.a);
  
  //     update shear and stratification
  meanflow_shear (&nl, &gotm_cnpar);
  meanflow_stratification (&nl, &gotm_buoy_method,
			   &dt, &gotm_cnpar,
			   turbulence_nuh.a, turbulence_gamh.a);

  //    compute turbulent mixing
  switch (turbulence_turb_method) {

  case 0:
    //        do convective adjustment
    meanflow_convectiveadjustment (&nl, turbulence_num.a, turbulence_nuh.a,
				   &turbulence_const_num, &turbulence_const_nuh,
				   &gotm_buoy_method, &meanflow_gravity,
				   &meanflow_rho_0);
    break;

  case 99: {
    //        update KPP model
    realtype precip = airsea_precip + airsea_evap;
    realtype tFlux, btFlux, sFlux, bsFlux;
    realtype tRad[nl+1], bRad[nl+1]; // fixme: check this
    util_convert_fluxes (&nl, &meanflow_gravity, &meanflow_cp,
			 &meanflow_rho_0,
			 &airsea_heat, &precip,
			 meanflow_rad.a, meanflow_t.a, meanflow_s.a,
			 &tFlux, &sFlux, &btFlux, &bsFlux,
			 tRad, bRad);

    kpp_do_kpp (&nl, &meanflow_depth,
		meanflow_h.a, meanflow_rho.a, meanflow_u.a,
		meanflow_v.a, meanflow_nn.a, meanflow_nnt.a, meanflow_nns.a,
		meanflow_ss.a,
		&meanflow_u_taus, &meanflow_u_taub,
		&tFlux, &btFlux, &sFlux, &bsFlux,
		tRad, bRad, &meanflow_cori);
    break;
  }
    
  default:
    //        update one-point models
    turbulence_do_turbulence (&nl, &dt, &meanflow_depth,
			      &meanflow_u_taus,
			      &meanflow_u_taub, &meanflow_z0s, &meanflow_z0b,
			      meanflow_h.a, meanflow_nn.a, meanflow_ss.a,
#ifdef SEAGRASS
			      xP.a
#else
			      NULL
#endif
			      );
  }
}

/**
The initialisation tries to call the relevant initialisation functions
from GOTM, as they appear in the original
[gotm.F90](https://github.com/gotm-model/code/blob/v5.2/src/gotm/gotm.F90#L129)
source code. */

event defaults (i = 0)
{
  // cnpar       [float, minimum = 0, maximum = 1]
  //             Constant for the theta scheme used for time integration of
  //             diffusion-reaction components. \theta=0.5 for Cranck-Nicholson
  //             (second-order accurate), \theta=0 for Forward Euler (first-
  //             order accurate), \theta=1 for Backward Euler (first-order
  //             accurate). Note that only \theta=1 guarantees positive
  //             solutions for positive definite systems.
  gotm_cnpar = 1.;

  // buoy_method [integer]
  //               method to compute mean buoyancy
  //               1: from equation of state (i.e. from potential temperature
  //                                          and salinity)
  //               2: from prognostic equation
  gotm_buoy_method = 1;

  meanflow_gravity = G;
  
  mtridiagonal_init_tridiagonal (&nl);
  realtype latitude = 0.;
  meanflow_post_init_meanflow (&nl, &latitude);
  if (turbulence_turb_method == 99) {
    integer namlst = -1;
    kpp_init_kpp (&namlst, "", &nl, &meanflow_depth, meanflow_h.a,
		  &G, &meanflow_rho_0);
  }
  else
    turbulence_post_init_turbulence (&nl);
}

/**
Vertical diffusion is part of the "viscous terms" of the multilayer
event loop. We just copy the relevant fields into the corresponding
GOTM fields, call GOTM and retrieve the updated fields. */

event viscous_term (i++)
{
  struct { realtype * x, * y; } gotm_u = { meanflow_u.a, meanflow_v.a };
  struct { realtype * x, * y; } gotm_t = { &airsea_tx, &airsea_ty };
  foreach() {
    foreach_dimension()
      *gotm_t.x = airsea_tau.x[];
    
    double z = zb[];
    foreach_layer() {
      meanflow_h.a[point.l + 1] = h[];
      foreach_dimension()
	gotm_u.x[point.l + 1] = u.x[];
      if (T.i >= 0)
	meanflow_t.a[point.l + 1] = T[];
      if (S.i >= 0)
	meanflow_s.a[point.l + 1] = S[];
      meanflow_z.a[point.l + 1] = z + h[]/2.;
      z += h[];
    }

    airsea_heat = airsea_heat_flux[];
    airsea_i_0 = airsea_swr_flux[];
    gotm_step (i);

    foreach_layer() {
      foreach_dimension()
	u.x[] = gotm_u.x[point.l + 1];
      if (T.i >= 0)
	T[] = meanflow_t.a[point.l + 1];
      if (S.i >= 0)
	S[] = meanflow_s.a[point.l + 1];
    }
  }
}

/**
The GOTM cleanup function deallocates Fortran fields etc. */

event cleanup (t = end)
{
  gotm_clean_up();
}

/**
## Utility functions

This function uses GOTM to initialize a vertical temperature profile
corresponding to a squared buoyancy frequency `NN` for a given constant
salinity and top temperature. */

void constant_NNT (double T_top, double S_const, double NN,
		   scalar T)
{
  foreach() {
    double z = zb[];
    foreach_layer()
      z += h[];
    T[0,0,nl-1] = T_top;
    for (int l = nl - 2; l >= 0; l--) {
      double pFace = G*(z - h[0,0,l]/2.);
      double alpha = eqstate_eos_alpha (&S_const, &T[0,0,l+1], &pFace, &G,
					&meanflow_rho_0);
      T[0,0,l] = T[0,0,l+1] - 1./(G*alpha)*NN*h[0,0,l];
      z -= h[0,0,l];
    }
  }
}
