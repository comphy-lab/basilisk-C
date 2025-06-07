/**
# Fingering instability in a colloidal solution of active microswimmers.
![64 _pullers_ in a periodic domain](fingering_without_penetration/fingering_without_penetration.mp4)(loop)

*/

/**
## Particle penetration treatment.
Provided that particle's have _zero_ radius, penetration could occur.
To avoid penetration, one classic idea would be to impose a body force
following a given potential (e.g. Lennar-Jones...).
This provides couple of problems in our configuration:
 - problem could become _stiff_ to solve;
 - I'm inevitably adding a body force to the swimmer.
Especially the last one, is intrinsically wrong.
My swimmer *must* be force free. Hence, adding a repulsive
potential force, I'm modifying the physics of the system
(and adding parasite forces...).

To avoid this, I'll try here a different approach: simply 
correcting the particle's location IF their distance (GAP)
is smaller than a given quantity.
It seems to work remarkably well.
It also avoids the _bounce-back_ effect that one would see
in an inertial flow. We are not inertial here.*/


/**
tHree
*/

#include "navier-stokes/centered.h"
#include "view.h"
/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */

double HEIGHT = 32.0;
double THETA = M_PI/2.*1;
double RADIUS = 1.0;
double B = 1.0; // reorientation owing to gravity

#include "fpicella/src/periodic-shift-treatment.h"

#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double theta; double thetan; double omega; double omega2; double r; \
coord B; double Thrust;

/** A pusher? A puller?
		draw your microswimmer! .*/
#define alpha +2
#define beta  +2

//#define Particle_Advection_ON 1
/** In this implementation, I will use a custom (and simpler...)
		timestepper for advecting the particle's location.
		In particular, simple Explicit Euler is used.
		It is way simpler than Antoon's Runge Kutta, but it does the job!*/
/** Most important, it couples well with the anti
		penetration alghoritm that I've employed.
		Here I just work on particle's location, no on velocity. Why?
		My suggestion (it's a model!) is that we're in Stokes regime.
		You can not _bounce_ when quasi-steady ;) */

/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */
//#define OMEGA_EQUATION() cos(p().theta)/(2.*B) + interpolate_linear(locate (x, y, z), omega, x, y, z)/2.  // FP 20250308 
#define OMEGA_EQUATION() 0.
#include "fpicella/src/tracer-particles-FP.h"
#include "fpicella/src/driver-tHree-smooth.h"

#include "Antoonvh/scatter2.h"


/**
	Vorticity field. */
scalar omega[];

/**
	Variable associated to particles that are modelling microswimmers */
Particles microswimmers;

FILE * output_file; // global output file pointer

scalar un[]; // reference velocity field, for stopping the simulation.

double NP = 8.; // number of agents I want to work with...


int main()
{
	
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (64. [0]);
  //DT = HUGE [0];
  DT = 1.0 [0];
  
  origin (-L0/2., -L0/2.);

  periodic (right);
  periodic (top);
  
  stokes = true;
  //TOLERANCE = 1e-3;

  //output_file = fopen ("Output.dat", "w");

  
//  for (N = 32; N <= 256; N *= 2)
//		for(RADIUS = 0.05; RADIUS <= 0.45; RADIUS += 0.05)
	N=64;
    	run();

	//fclose(output_file);
}

#define WIDTH 0.5
#define EPS 1e-14

event init (t = 0) {
/**
	Initialize microswimmer particles. 
	When using init_tp_square, the total number
	of particles is NP*NP.*/
	microswimmers = init_tp_square(NP);
	foreach_particle_in(microswimmers){
		p().B.y=0.*-1.;
		p().r = RADIUS;
		p().Thrust = 4.;
		p().theta = THETA+0.*noise()*M_PI*2.; 
	}

  /**
  Viscosity is  unity. */
  
  mu = fm;

  ///**
  //The channel geometry is defined using Constructive Solid Geometry. */  
  //mask (y > +HEIGHT/2. ? top : none);
  //mask (y < -HEIGHT/2. ? bottom : none);

  //mask (x > +HEIGHT/2. ? right : none);
  //mask (x < -HEIGHT/2. ? left : none);
	
	u.n[top]    = dirichlet(0.);
	u.t[top]    = dirichlet(0.);
	u.n[bottom] = dirichlet(0.);
	u.t[bottom] = dirichlet(0.);
	u.t[right]  = dirichlet(0.); 
	u.n[right]  = dirichlet(0.);
	u.n[left]   = dirichlet(0.);
	u.t[left]   = dirichlet(0.);
	
	//* Initialize reference velocity field.*/
	  foreach()
    un[] = u.y[];
}


event log_make_video (t += 1.0) {
  squares ("omega", linear = true);
  scatter (microswimmers, pc = {256, 256, 256});
  isoline ("bodyPlot", 0.25, lc = {256,256,256}, lw = 2);
  vectors (u = "u", scale = 0.5, lc = {256,256,256}, level = 8);
  box();
  save ("fingering_without_penetration.mp4");
}

/** End event, to properly unallocate variables.*/
event end (t = 100);

event acceleration (i++){
	compute_microswimmer_forcing_smooth(microswimmers);
	//compute_repulsion(microswimmers,4.0,false,false,false,false,6.); // look at the driver-tHree file
																																		// to know what argument is meant
																																		// to do ;)
	a = av; // assign acceleration term. That's capital!
}

event properties(i++){
	vorticity(u,omega);
	compute_bodyPlot(microswimmers);
}

/** Microswimmer's kinematics, using simple Forward Euler.*/
event end_timestep(i++){
	foreach_particle(){
	/** Standard behaviour.*/
		foreach_dimension()
			p().x += p().u.x*dt; // The most simple timestepper possible:
													 // Explicit Euler.
  	/** Detect contact (GAP<0) and correct kinematics..*/
  	for (int _k_particle = 0; _k_particle < pn[_l_particle]; _k_particle++) {
  	  if(_k_particle != _j_particle){
  	    coord DISTANCE; // relative distance
  	    foreach_dimension()
  	    DISTANCE.x = pl[_l_particle][_k_particle].x - pl[_l_particle][_j_particle].x;
  	    double distance = sqrt(sq(DISTANCE.x)+sq(DISTANCE.y));
        double angle = atan2(DISTANCE.y/distance,DISTANCE.x/distance);
				double GAP = distance - 6.0*p().r;
				/** What if the GAP < 0?
				p() are getting too close, action required!*/
			  if(GAP<0.){
					/** Correct position to avoid penetration.*/
					/** Correction is simply proportional to the amount
							of penetration (GAP<0) that is observed.*/
					p().x += GAP*cos(angle);
					p().y += GAP*sin(angle);
					/** Since penetration is symmetric between the two
							particles, the new outcome GAP will be symmetrical
							distributed to the two particles as well.*/
			  }
  	  }
  	}
	}
}