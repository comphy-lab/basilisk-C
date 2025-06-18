/**
tHree
*/

#include "navier-stokes/centered.h"
#include "view.h"
/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */

double NP = 2.; // number of agents I want to work with...
double HEIGHT =128.0;
double THETA = M_PI/2.*1;
double RADIUS = 1.0;
double B = 10.0; // reorientation owing to gravity
double STERIC_GAP = 5.0; // interaction GAP, before an action is taken,
												 // to avoid steric interaction.
												 //
												 //
double BINARY_DISTANCE = 10.;
double BUOYANCY = 0.;

#include "fpicella/src/periodic-shift-treatment.h"

#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double theta; double thetan; double omega; double omega2; double r; \
coord B; double Thrust;

/** A pusher? A puller?
		draw your microswimmer! .*/
#define alpha +2.0
#define beta  +2.0

/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */
//#define OMEGA_EQUATION() cos(p().theta)/(2.*B) + interpolate_linear(locate (x, y, z), omega, x, y, z)/2.  // FP 20250308 
#define OMEGA_EQUATION() interpolate_linear(locate (x, y, z), omega, x, y, z)/2.  // FP 20250308 
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

  periodic (right);
  //periodic (top);
  
  stokes = true;
  //TOLERANCE = 1e-3;

  output_file = fopen ("Output.dat", "w");

  
//  for (N = 32; N <= 256; N *= 2)
//		for(RADIUS = 0.05; RADIUS <= 0.45; RADIUS += 0.05)
	N=128;
	BINARY_DISTANCE = L0/2.;
		for(BINARY_DISTANCE = L0/2.; BINARY_DISTANCE>=RADIUS; BINARY_DISTANCE*=0.9)
    	run();

	fclose(output_file);
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
		p().Thrust = 4.;
		p().theta = THETA+noise()*M_PI*0.; 
/** Random initialization .*/
		p().x = -BINARY_DISTANCE/2. +_j_particle*p().r*BINARY_DISTANCE;
		p().y = 0.;
	}

  /**
  Viscosity is  unity. */
  
  mu = fm;

  ///**
  //The channel geometry is defined using Constructive Solid Geometry. */  
  mask (y > +HEIGHT/2. ? top : none);
  mask (y < -HEIGHT/2. ? bottom : none);

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


//event log_make_video (i += 1) {
//  squares ("omega", linear = true);
//  scatter (microswimmers, pc = {256, 256, 256});
//  isoline ("bodyPlot", 0.25, lc = {256,256,256}, lw = 2);
//  vectors (u = "u", scale = 0.25, lc = {256,256,256}, level = 8);
//  box();
//  save ("movie.mp4");
//}

event acceleration (i++){
	compute_microswimmer_forcing_smooth(microswimmers);
	a = av; // assign acceleration term. That's capital!
}

event properties(i++){
	vorticity(u,omega);
	compute_bodyPlot(microswimmers);
}

/** Microswimmer's kinematics, is turned off.*/

/** When computation has conveged, write the output file.*/
event logfile ( i += 1; i <= 1000) {
  double du = change (u.y, un);
  if (i > 5 && du < 1e-5){ // since I'm looking for a steady solution, this must be quite small...
		foreach_particle_in(microswimmers){
  		fprintf(output_file,"%+3.2e %+04i %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+04i\n",
                L0, N, RADIUS, alpha, beta, p().Thrust, BUOYANCY, B, BINARY_DISTANCE, p().u.x,p().u.y,OMEGA_EQUATION(),_j_particle);
  		fflush(output_file);
		}
    return 1; /* stop */
  }
}

/**
# Binary interactions between tHree microswimmers
~~~gnuplot Binary interactions between two tHree microswimmers
set grid
set xlabel "d/R"
set title "binary interactions, spanwise aligned, theta = pi/2. (poiting y direction)"

# Plot data from file
plot 'Output.dat' using ( $13==0 && $7==0 ? ($9/$3) : 1/0 ):( $13==0  && $7==0 ? +$11 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "green" title "u_y",\
     'Output.dat' using ( $13==0 && $7==0 ? ($9/$3) : 1/0 ):( $13==0  && $7==0 ? +$10 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "blue" title "u_x ",\
     'Output.dat' using ( $13==0 && $7==0 ? ($9/$3) : 1/0 ):( $13==0  && $7==0 ? +$12 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "red" title "omega"
~~~
*/
