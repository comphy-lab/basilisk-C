double nu     = 1.; // fluid viscosity;
double nuRatio= 100.0;// particle-fluid viscosity ratio
double deltaT = 0.1;
double PROPULSION = 0.;
double SEDIMENTATION = 1.;
double RADIUS = 0.25;
double HEIGHT = 1.0+1e-2; // channel height, for the confinement
#define channel difference(-x+HEIGHT/2.,-x-HEIGHT/2.)

/*
Parameters to play with particles...
*/
#define NPARTICLES 1 
/*
Features of the forces...
*/
#define ALPHA 1.5 // sideways asymmetry
#define BETA  1.5 // front (+) back (-) asymmetry
//#define PROPULSION 5.0 // strenth of PROPULSION
//#define SEDIMENTATION 0.0 // strength of SEDIMENTATION
#define B 1e10 // gravitaxis orientation // smaller value? Stronger reorientation!

// If embedded boundaries are needed (fancy confinement...)
#include "ghigo/src/myembed.h"
#include "fpicella/src/centered_embed_variable_viscosity.h"

//// If embedded boundaries are not required
//scalar cs[];
//face vector fs[];
//#include "fractions.h"
//#include "navier-stokes/centered.h"

/**
## Setup
We need a field for viscosity so that the embedded boundary metric can
be taken into account. */
face vector muv[]; // variable viscosity
face vector av[] ; // varialbe acceleration...accounting for microswimmer's forces
#include "fpicella/src/periodic-shift-treatment.h"


//#define Particle_Advection_ON 1
#include "fpicella/src/microswimmer_particle_driver.h"

#include "view.h"

/**
We define the mesh adaptation parameters. */

int lmin = 5;  // Min mesh refinement level (l=6 is 2pt/d)
int lmax = 5; // Max mesh refinement level (l=10 is 36pt/d, l=13 is 292pt/d)
#define cmax (1.e-6) // Absolute refinement criteria for the velocity field

int main()
{
	display_control(deltaT,0.001,100);
	display_control(nuRatio,1,1000);
	display_control(RADIUS,0.1,0.45);
	display_control(SEDIMENTATION,-100,100);
	display_control(PROPULSION,-100,100);
	display_control(lmax,lmin,10);
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (2 [0]);
  DT = 1.0 [0];

  origin (-L0/2., -L0/2.);

  stokes = true;

  periodic(left);
  periodic(top);

  N = 1 << (lmax);//initialize at maximum refinement, so to have _child_

  init_grid (N);

  run();  

}

event init(i=0)
{
	u.t[left]   = dirichlet(0.);
	u.n[left]   = dirichlet(0.);
	u.t[right]  = dirichlet(0.);
	u.n[right]  = dirichlet(0.);
	u.t[top]   = dirichlet(0.);
	u.n[top]   = dirichlet(0.);
	u.t[bottom]  = dirichlet(0.);
	u.n[bottom]  = dirichlet(0.);

#if EMBED
#if TREE

  astats ss;
  int ic = 0;
  do {
    ic++;
		solid(cs,fs,channel);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
                        maxlevel = (lmax), minlevel = (lmin));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
	u.t[embed] = dirichlet(0.);
	u.n[embed] = dirichlet(0.);
	uf.t[embed] = dirichlet(0.);
	uf.n[embed] = dirichlet(0.);
#endif

  foreach_face()
    muv.x[] = (nu)*fm.x[];
  boundary ((scalar *) {muv});

//  for (scalar s in {u, p, pf})
//    s.third = true;
}

event logfile (t += deltaT; i <= 10000) {
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

event properties(i++){
	foreach_particle(){
		p().r = RADIUS;
		p().x = 0.;
		p().y = 0.;
	}
}

event adapt (i++) // and not classic adapt, so to keep mesh identical during subiterations...
{
	solid(cs,fs,channel);
  adapt_wavelet ({u}, (double[]) {(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (lmin));
}

