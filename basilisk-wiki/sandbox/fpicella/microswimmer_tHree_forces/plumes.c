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
//#define OMEGA_EQUATION() 0.
//#define OMEGA_EQUATION() interpolate_linear (locate (x, y, z),omega, x, y, z)

double HEIGHT = 32.0;
double THETA = M_PI/2.*1;
double RADIUS = 1.0;
double FACTOR = 1.;
double B = 10.; // reorientation owing to gravity

#include "fpicella/src/periodic-shift-treatment.h"

#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double theta; double thetan; double omega; double omega2; double r; \
coord B; double Thrust;

#define alpha +2
#define beta  +2

#define Particle_Advection_ON 1
#define OMEGA_EQUATION() cos(p().theta)/(2.*B) + interpolate_linear(locate (x, y, z), omega, x, y, z)/2.  // FP 20250308 
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

double NP = 7.;


int main()
{
	
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (64. [0]);
  DT = HUGE [0];
  
  origin (-L0/2., -L0/2.);

  periodic (right);
  //periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-3;

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
	Initialize microswimmer particles. */
	//microswimmers = init_tp_circle(NP);
	microswimmers = init_tp_square(5);
	foreach_particle_in(microswimmers){
		p().B.y=-1.;
		p().r = RADIUS;
		p().Thrust = 4.;
		//p().y = alternating_series(_j_particle)*9.*p().r;
		//p().x = 0.;
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


event logfile (t += 0.125; i <= 5000) {
  squares ("omega", linear = true);
  scatter (microswimmers, pc = {256, 256, 256});
  isoline ("bodyPlot", 0.25, lc = {256,256,256}, lw = 2);
  vectors (u = "u", scale = 0.5, lc = {256,256,256}, level = 8);
  box();
  save ("plumes.mp4");
}

event acceleration (i++){
	compute_microswimmer_forcing_smooth(microswimmers);
	compute_repulsion(microswimmers,4.0,true,true,true,true,6.);
	a = av;
}

event properties(i++){
	vorticity(u,omega);
	//foreach_particle_in(microswimmers)
	//	p().theta = THETA;
/** For plotting purposes only...*/
	compute_bodyPlot(microswimmers);
}

/**
# Plumes of bioconvection (?)

![Some plumes from heart-shaped microswimmers](plumes/plumes.mp4)

*/