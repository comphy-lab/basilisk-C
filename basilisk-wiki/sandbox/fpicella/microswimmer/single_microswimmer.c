double nu     = 1.; // fluid viscosity;
double nuRatio= 100.0;// particle-fluid viscosity ratio
double deltaT = 1.0;

/*
Parameters to play with particles...
*/
#define RADIUS 1.0
#define NPARTICLES 2 
/*
Features of the forces...
*/
#define ALPHA 1.5 // sideways asymmetry
#define BETA  1.5 // front (+) back (-) asymmetry
#define PROPULSION 1.0 // strenth of PROPULSION
#define SEDIMENTATION 0.0 // strength of SEDIMENTATION

// If embedded boundaries are needed (fancy confinement...)
//#include "ghigo/src/myembed.h"
//#include "fpicella/src/centered_embed_variable_viscosity.h"

// If embedded boundaries are not required
#include "fractions.h"
#include "navier-stokes/centered.h"

/**
## Setup
We need a field for viscosity so that the embedded boundary metric can
be taken into account. */
face vector muv[]; // variable viscosity
face vector av[] ; // varialbe acceleration...accounting for microswimmer's forces
#include "fpicella/src/periodic-shift-treatment.h"

#include "fpicella/src/microswimmer_particle_driver.h"

#include "view.h"

/**
We define the mesh adaptation parameters. */

int lmin = 4;  // Min mesh refinement level (l=6 is 2pt/d)
int lmax = 6; // Max mesh refinement level (l=10 is 36pt/d, l=13 is 292pt/d)
#define cmax (1.e-2) // Absolute refinement criteria for the velocity field

int main()
{
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (16 [0]);
  DT = 1.0 [0];

  origin (-L0/2., -L0/2.);

  stokes = true;

  periodic(left);
  //periodic(top);

  N = 1 << (lmax);//initialize at maximum refinement, so to have _child_

  init_grid (N);

  run();  

}

event init(i=0)
{
	u.t[top]     = dirichlet(+1.0);
	u.n[top]     = dirichlet(0.);
	u.t[bottom]  = dirichlet(-1.0);
	u.n[bottom]  = dirichlet(0.);

  foreach_face()
    muv.x[] = (nu)*fm.x[];
  boundary ((scalar *) {muv});

//  for (scalar s in {u, p, pf})
//    s.third = true;
}

event logfile (t += deltaT; i <= 1000) {
	foreach_particle()
		fprintf(stderr,"%+6.5e %+6.5e %+6.5e %+6.5e \n",p().x,p().y,p().theta,p().areaCenter);
}

event properties(i++){
	/*
	muv face field is calculated so to "stiffen" the fluid at the location of each microswimmer
	computation is done within microswimmer_particle_driver.h
	*/
	event("compute_variable_viscosity_field");
  mu = muv;
  boundary ((scalar *) {muv});
	/*
	Compute the forces required to model the microswimmer motion
	*/
	event("compute_propulsion_forces");
	/*
	Compute vorticity field...to update particle's velocity and locations
	*/
	vorticity(u,omega);
}

