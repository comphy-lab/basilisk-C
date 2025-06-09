/**
# Single tHree micro-swimmer

![Single heart-shaped microswimmers in a periodic domain](single_swimmer_base_case/single_swimmer.mp4)(loop)
Of course, this aims to mimick a [*puller*](https://doi.org/10.1146/annurev-fluid-121021-042929)...

*/



int alternating_series(int n) {
    if (n == 0) return 0;
    int k = (n + 1) / 2; // 1, 1, 2, 2, 3, 3, ...
    return (n % 2 == 1) ? k : -k;
}

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
double FACTOR = 1.;
double B = 1.; // reorientation owing to gravity

#include "fpicella/src/periodic-shift-treatment.h"

#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double theta; double thetan; double omega; double omega2; double r; \
coord B; double Thrust;

#define alpha +2
#define beta  +2

#define Particle_Advection_ON 1
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

double NP = 1.; // number of agents I want to work with...


int main()
{
	
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (32. [0]);
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
	N=32;
    	run();

	//fclose(output_file);
}

#define WIDTH 0.5
#define EPS 1e-14

event init (t = 0) {
/**
	Initialize microswimmer particles. */
	//microswimmers = init_tp_circle(NP);
	microswimmers = init_tp_square(NP);
	foreach_particle_in(microswimmers){
		p().B.y=0.*-1.;
		p().r = RADIUS;
		p().Thrust = 4.;
		//p().x = L0/2.-p().r*4.+alternating_series(_j_particle)*8.*p().r;
		//p().x = L0/2.+_j_particle*L0/NP*p().r;
		//double SPACING = 8.;
		//p().x = -L0/2.+(-SPACING/2.+_j_particle*SPACING)*p().r;
		//p().y = alternating_series(_j_particle)*9.*p().r;
		//p().y = 0.;
		// forces will be normalized by the cylinder area
		// this is done in driver-tHree-smooth.h
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
  save ("single_swimmer.mp4");
}

/** End event, to properly unallocate variables.*/
event end (t = 100);

event acceleration (i++){
	compute_microswimmer_forcing_smooth(microswimmers);
	compute_repulsion(microswimmers,3.0,false,false,false,false,12.); // look at the driver-tHree file
																																		// to know what argument is meant
																																		// to do ;)
	a = av; // assign acceleration term. That's capital!
}

event properties(i++){
	vorticity(u,omega);
	compute_bodyPlot(microswimmers);
}

/**
# Velocity for a single tHree swimmer.
~~~gnuplot tHree swimmer
set grid
set xlabel "t"
set ylabel "U_y"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot 'particle_000.dat' using ($1):($6) with points pt 7 ps 2 lc rgb "green" title "present"
~~~
I should observe:

- instantaneously getting to asymptotic value
- no influence of the particle position on it's velocity
- no jumps when crossing a boundary (I'm periodic...)
- I should observe no drift in time, since I've got a force-free system.

Still, some oscillations due to lagrangian-eulerian interpolation.
but honestly, it's a piece of cake.
*/
