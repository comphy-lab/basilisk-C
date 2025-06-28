/**
# Swimming velocity of a single, non force-free (buoyant) particle.
Placed between two plates.
(like an infinite vertical pipe...)
![Single heart-shaped microswimmers in a y periodic domain + BUOYANCY (non-force-free!)](single_y_periodic_buoyancy/single_y_periodic_buoyancy.mp4)(loop)

### Non-periodic boundaries
are necessary when buoyancy is on.
Why?
Because the microswimmer + buoyancy IS NOT FORCE FREE.
Hence, I would observe some _run-out_...
Simple no-slip walls will do the job, and guarantee a 
zero velocity _far away_ from the swimmer.
*/

#include "navier-stokes/centered.h"
#include "view.h"
/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */

double NP = 1.; // number of agents I want to work with...
double HEIGHT =64.0;
double THETA = M_PI/2.*1;
double RADIUS = 1.0;
double THRUST = 4.;
double BUOYANCY = -4*0.05;
double B = 1.0; // reorientation owing to gravity
double STERIC_GAP = 4.0; // interaction GAP, before an action is taken,
												 // to avoid steric interaction.

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
#define OMEGA_EQUATION() cos(p().theta)/(2.*B) + interpolate_linear(locate (x, y, z), omega, x, y, z)/2.  // FP 20250308 
//#define OMEGA_EQUATION() 0.
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

int main()
{
	
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (64. [0]);
  //DT = HUGE [0];
  DT = 1.0 [0];
  
  origin (-L0/2., -L0/2.);

  //periodic (right);
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
  Boundary condition on top-bottom boundaries.*/
	u.n[top] = dirichlet(0.);
	u.t[top] = dirichlet(0.);
	u.n[bottom] = dirichlet(0.);
	u.t[bottom] = dirichlet(0.);

/**
	Initialize microswimmer particles. 
	When using init_tp_square, the total number
	of particles is NP*NP.*/
	//microswimmers = init_tp_square(NP);
	microswimmers = init_tp_circle(NP);
	foreach_particle_in(microswimmers){
/** This is where I set a non-zero buoyancy force.*/
		p().B.y=BUOYANCY;
		p().r = RADIUS;
		p().Thrust = THRUST;
		p().theta = THETA;//+noise()*M_PI*2.; 
/** Random initialization .*/
		p().x =     L0/2.*noise();
		p().y = HEIGHT/2.*noise();
		if(_j_particle==0){
			p().x = 0.;
			p().y = 0.;
		}
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


event log_make_video (i += 10) {
  squares ("omega", linear = true);
  scatter (microswimmers, pc = {256, 256, 256});
  isoline ("bodyPlot", 0.25, lc = {256,256,256}, lw = 2);
  vectors (u = "u", scale = 0.25, lc = {256,256,256}, level = 8);
  box();
  save ("single_y_periodic_buoyancy.mp4");
}

/** End event, to properly unallocate variables.*/
event end (t = 2000);

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
	/** Angular kinematics.*/
		p().omega = OMEGA_EQUATION();
		p().theta += p().omega*dt;

//		/** Detect if particle is in contact with top or bottom boundary.*/
//		/** Hence, that's simply a matter of y-aligned dynamics.*/
//		if(p().y+p().r*STERIC_GAP*2.>HEIGHT/2.) // impacting with top boundary
//			p().y = HEIGHT/2.-p().r*STERIC_GAP*2.;
//		if(p().y-p().r*STERIC_GAP<-HEIGHT/2.) // impacting with bottom boundary
//			p().y = -HEIGHT/2.+p().r*STERIC_GAP;
//
		/** Detect if particle is in contact with top or bottom boundary.*/
		/** Hence, that's simply a matter of y-aligned dynamics.*/
		/** IN THE X DIRECTION NOW! */
		if(p().x+p().r*STERIC_GAP*2.>HEIGHT/2.) // impacting with top boundary
			p().x = HEIGHT/2.-p().r*STERIC_GAP*2.;
		if(p().x-p().r*STERIC_GAP<-HEIGHT/2.) // impacting with bottom boundary
			p().x = -HEIGHT/2.+p().r*STERIC_GAP;



  	/** Detect particle-particle contact (GAP<0) and correct kinematics..*/
  	for (int _k_particle = 0; _k_particle < pn[_l_particle]; _k_particle++) {
  	  if(_k_particle != _j_particle){
  	    coord DISTANCE; // relative distance
  	    foreach_dimension()
  	    DISTANCE.x = pl[_l_particle][_k_particle].x - pl[_l_particle][_j_particle].x;
  	    double distance = sqrt(sq(DISTANCE.x)+sq(DISTANCE.y));
        double angle = atan2(DISTANCE.y/distance,DISTANCE.x/distance);
				double GAP = distance - STERIC_GAP*p().r;
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

/**
~~~gnuplot Swimming velocity of a negatively buoyant microswimmer
reset
set title "Swimming velocity, active microswimmer + buoyancy."
set xlabel "time"
set ylabel "U_{swim}"
set grid
#set term wxt enhanced persist

# Set the moving average window size here
window = 50

# Correct: no backslashes for $1/$6
avg_cmd = sprintf("<awk -v w=%d '{x[NR]=$1; y[NR]=$6} NR>=w {sum=0; for(i=NR-w+1;i<=NR;i++) sum+=y[i]; print x[NR], sum/w}' particle_000.dat", window)

plot "particle_000.dat" using 1:6 with points title "Raw Data", \
     avg_cmd using 1:2 with lines linewidth 4 title sprintf("Moving Avg (%d)", window)
~~~
*/
