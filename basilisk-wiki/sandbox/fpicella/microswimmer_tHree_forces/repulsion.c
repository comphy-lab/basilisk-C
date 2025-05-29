/**
# Single cylindrical particle between flat plates.
Velocity as a function of confinement, for constant force.
No use of embed. _Complex_ boundary is obtained using *maks*.
*/

#include "navier-stokes/centered.h"
#include "view.h"
/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */
#define OMEGA_EQUATION() 0.
//#define Particle_Advection_ON 1

double HEIGHT = 1.0;
double THETA = M_PI/2.*3;
double RADIUS = 0.25;
double FACTOR = 1.;

#include "fpicella/src/periodic-shift-treatment.h"

#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double theta; double thetan; double omega; double omega2; double r; \
coord B; double Thrust;

#define alpha +2
#define beta  +2

#define Particle_Advection_ON 1
#include "fpicella/src/tracer-particles-FP.h"
#include "fpicella/src/driver-tHree-smooth.h"


/**
	Vorticity field. */
scalar omega[];

/**
	Variable associated to particles that are modelling microswimmers */
Particles microswimmers;

FILE * output_file; // global output file pointer

scalar un[]; // reference velocity field, for stopping the simulation.

double NP = 1.;


int main()
{
	display_control(THETA,-10,10);
	
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (16. [0]);
  DT = HUGE [0];
  
  origin (-L0/2., -L0/2.);

  //periodic (right);
  //periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-7;

  output_file = fopen ("Output.dat", "w");

	RADIUS = 0.25;
  
//  for (N = 32; N <= 256; N *= 2)
//		for(RADIUS = 0.05; RADIUS <= 0.45; RADIUS += 0.05)
	N=64;
    	run();

	fclose(output_file);
}

#define WIDTH 0.5
#define EPS 1e-14

event init (t = 0) {
/**
	Initialize microswimmer particles. */
	microswimmers = init_tp_circle(NP);
	foreach_particle(){
		//p().x = 0.;
		//p().y = 0.;
		p().B.y=-1.;
		p().r = RADIUS;
		p().Thrust = 0.;
		// forces will be normalized by the cylinder area
		// this is done in driver-tHree-smooth.h
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


/**
We check for a stationary solution. */

event logfile (t += 0.05; i <= 2000) {
//  double du = change (u.y, un);
	foreach_particle_in(microswimmers){
		fprintf(stderr,"PARTICLE %+6.5e %+6.5e %+6.5e %+6.5e \n",p().x,p().y,p().u.x,p().u.y);
  	//if (i > 10 && du < 1e-6){
		//	fprintf(output_file,"%+04d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e \n",N,p().r,p().x,p().y,p().u.x,p().u.y,HEIGHT);
		//	fflush(output_file);
    //	return 1; /* stop */
		//}
	}
}

event acceleration (i++){
	compute_microswimmer_forcing_smooth(microswimmers);
	compute_repulsion(microswimmers);
	a = av;
}

event properties(i++){
	foreach_particle_in(microswimmers)
		p().theta = THETA;
}
