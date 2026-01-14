/**
# Constant properties
This file contains functions to update physical properties in simulations 
of two-phase flow with phase change in porous media. We want to solve the
Navier-Stokes equations using constant properties for the gas phase throughout
the domain, both in the gas and in the porous solid regions.
To achieve this, we overload the `properties` event to set the density and viscosity
of the gas also in the solid regions. Additionally, we set the thermal conductivity
and mass diffusivity for multicomponent simulations.
*/

#ifndef F_ERR
# define F_ERR 1.e-10
#endif

#define pavg(p,vg,vs) (p*vg+(1-p)*vs)

extern double rhoG, rhoS;
extern double muG;
extern face vector alphav;
extern scalar rhov;
extern scalar f;
extern scalar porosity;

/**
# Thermal conductivity boundary conditions
Given that we also want to account for different thermal conductivities
in different directions, lambda1v and lambda2v are defined as vectors.
This set of boundary conditions are needed to ensure the correct calculation
of the heat flux at the boundaries.
*/

#ifdef SOLVE_TEMPERATURE
extern vector lambda1v, lambda2v;
lambda2v.n[bottom] = neumann(0.); 
lambda2v.t[bottom] = neumann(0.);
lambda2v.n[left] = neumann(0.);
lambda2v.t[left] = neumann(0.);
lambda2v.n[right] = neumann(0.);
lambda2v.t[right] = neumann(0.);
lambda2v.n[top] = neumann(0.);
lambda2v.t[top] = neumann(0.);

lambda1v.n[bottom] = neumann(0.);
lambda1v.t[bottom] = neumann(0.);
lambda1v.n[left] = neumann(0.);
lambda1v.t[left] = neumann(0.);
lambda1v.n[right] = neumann(0.);
lambda1v.t[right] = neumann(0.);
lambda1v.n[top] = neumann(0.);
lambda1v.t[top] = neumann(0.);

extern double lambdaG, lambdaS;
extern double cpG, cpS;
#endif
#ifdef MULTICOMPONENT
extern unsigned int NGS;
extern scalar* DmixGList_S;
extern scalar* DmixGList_G;
#endif

/**
# Update properties with constant gas properties
This function updates the physical properties used in the Navier-Stokes equations
to have constant gas properties throughout the domain.
*/
void update_properties_constant (void) {
  scalar mu_centered[];
  foreach()
    mu_centered[] = muG/(porosity[]*f[] + (1. - f[]));

  foreach_face() {
    alphav.x[] = fm.x[]/rhoG;
    {
      face vector muv = mu;
      muv.x[] = face_value(mu_centered, 0)*fm.x[];
    }
  }

  foreach()
    rhov[] = rhoG*cm[];

#ifdef SOLVE_TEMPERATURE
  /**
  The pseudophase thermal conductivities are calcuated as weighted averages
  based on the porosity field.
  */
  foreach() {
    foreach_dimension() {
      lambda1v.x[] = f[] > F_ERR ? pavg (porosity[]/f[], lambdaG, lambdaS) : 0.;
      lambda2v.x[] = f[] < 1.-F_ERR ? lambdaG : 0.;
    }
  }
#endif

#ifdef MULTICOMPONENT
  double Dmixv =  2.05e-5; //Diff of CO in N2 at 500K, 1 atm

  foreach() {
    // set the same for all species
    for (int jj = 0; jj < NGS; jj++) {
      scalar Dmix2 = DmixGList_G[jj];
      Dmix2[] = f[] < 1. - F_ERR ? Dmixv : 0.;
    }
    for (int jj = 0; jj < NGS; jj++) {
      scalar Dmix2 = DmixGList_S[jj];
      /**
      In the solid region, the mass diffusivity is modified to account for
      the effect of the solid matrix.
      The simplest model is used here, where the diffusivity is scaled by the porosity.
      */
      Dmix2[] = f[] > F_ERR ? Dmixv * pow(porosity[], 4./3.) : 0.;
    }
  }
#endif
}

/**
# Overload properties event
Every time the properties event is called, we update the properties
to our desired constant gas properties.
*/

event properties (i++) {
  update_properties_constant();
}
