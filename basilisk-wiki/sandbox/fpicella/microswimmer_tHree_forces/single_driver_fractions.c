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
#define Particle_Advection_ON 1

double HEIGHT = 1.0;
double STRENGTH = 100.0;

#include "fpicella/src/periodic-shift-treatment.h"
#include "fpicella/src/tracer-particles-FP.h"
#include "fpicella/src/driver-tHree.h"


/**
	Vorticity field. */
scalar omega[];

/**
	Variable associated to particles that are modelling microswimmers */
Particles microswimmers;

double NP = 1.;

int main()
{
	display_control(STRENGTH,-1000,1000);
	
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (4. [0]);
  DT = HUGE [0];
  
  origin (-L0/2., -L0/2.);

  periodic (right);
  periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-7;
  
  //for (N = 16; N <= 64; N *= 2)
	N=64;
    run();
}

#define WIDTH 0.5
#define EPS 1e-14

event init (t = 0) {
/**
	Initialize microswimmer particles. */
	microswimmers = init_tp_circle(NP);

  /**
  Viscosity is  unity. */
  
  mu = fm;

  ///**
  //The channel geometry is defined using Constructive Solid Geometry. */  
  //mask (y > +HEIGHT/2. ? top : none);
  //mask (y < -HEIGHT/2. ? bottom : none);
	
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

event logfile (t += 0.05; i <= 1000) {
	foreach_particle_in(microswimmers){
		fprintf(stderr,"PARTICLE %+6.5e %+6.5e %+6.5e %+6.5e \n",p().x,p().y,p().u.x,p().u.y);
	}
}

event acceleration (i++){
	compute_microswimmer_forcing(microswimmers);
	a = av;
}
