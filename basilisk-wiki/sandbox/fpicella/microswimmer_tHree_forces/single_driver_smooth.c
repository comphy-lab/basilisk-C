/**
# Poiseuille flow in a periodic channel
but oriented in x direction.
I try to work with simple maks function
*/

#include "navier-stokes/centered.h"
#include "view.h"
/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */
#define OMEGA_EQUATION() 0.
//#define Particle_Advection_ON 1

double HEIGHT = 1.0;
double THETA = M_PI/2.*2;
double RADIUS = 0.25;
double FACTOR = 1.;

#include "fpicella/src/periodic-shift-treatment.h"

#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double theta; double thetan; double omega; double omega2; double r; \
coord B; double Thrust;

#define alpha +2
#define beta  +2

#include "fpicella/src/tracer-particles-FP.h"
#include "fpicella/src/driver-tHree-smooth.h"


/**
	Vorticity field. */
scalar omega[];

/**
	Variable associated to particles that are modelling microswimmers */
Particles microswimmers;

double NP = 1.;

int main()
{
	display_control(THETA,-10,10);
	
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (4. [0]);
  DT = HUGE [0];
  
  origin (-L0/2., -L0/2.);

  //periodic (right);
  //periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-7;
  
  //for (N = 16; N <= 64; N *= 2)
	N=32;
    run();
}

#define WIDTH 0.5
#define EPS 1e-14

event init (t = 0) {
/**
	Initialize microswimmer particles. */
	microswimmers = init_tp_circle(NP);
	foreach_particle(){
		p().x = 0.;
		p().y = 0.;
		p().B.x=0.;
		p().r = RADIUS;
		p().Thrust = 100.;
	}

  /**
  Viscosity is  unity. */
  
  mu = fm;

  ///**
  //The channel geometry is defined using Constructive Solid Geometry. */  
  //mask (y > +HEIGHT/2. ? top : none);
  //mask (y < -HEIGHT/2. ? bottom : none);

//  mask (x > +HEIGHT/2. ? right : none);
//  mask (x < -HEIGHT/2. ? left : none);
	
	u.n[top]    = dirichlet(0.);
	u.t[top]    = dirichlet(0.);
	u.n[bottom] = dirichlet(0.);
	u.t[bottom] = dirichlet(0.);
	u.t[right]  = dirichlet(0.); 
	u.n[right]  = dirichlet(0.);
	u.n[left]   = dirichlet(0.);
	u.t[left]   = dirichlet(0.);
}


/**
We check for a stationary solution. */

event logfile (t += 0.05; i <= 2000) {
	foreach_particle_in(microswimmers){
		fprintf(stderr,"PARTICLE %+6.5e %+6.5e %+6.5e %+6.5e \n",p().x,p().y,p().u.x,p().u.y);
	}
}

event acceleration (i++){
	compute_microswimmer_forcing_smooth(microswimmers);
	a = av;
}

event properties(i++){
	foreach_particle_in(microswimmers)
		p().theta = THETA;
}
