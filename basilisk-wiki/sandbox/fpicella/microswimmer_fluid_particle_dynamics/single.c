/**
# Poiseuille flow in a periodic channel
but oriented in x direction.
I try to work with simple maks function
*/

#include "navier-stokes/centered.h"
#include "view.h"
/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */
//#define OMEGA_EQUATION() 0.
#define OMEGA_EQUATION() interpolate_linear(locate (x, y, z),omega, x, y, z)
/**
Uncomment if you want particles to be advected */
#define Particle_Advection_ON 1

double HEIGHT = 1.0;
double STRENGTH = 100.0;
double ETA = 100.;
double RADIUS = 0.25;

#include "fpicella/src/periodic-shift-treatment.h"
#include "fpicella/src/tracer-particles-FP.h"
#include "fpicella/src/driver-FluidParticleDynamics.h"

/**
	Variable viscosity. */
face vector muv[];
/**
	Variable acceleration. */


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
	display_control(ETA,-1000,1000);
	display_control(RADIUS,-1000,1000);
	
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (4. [0]);
  DT = HUGE [0];
  
  origin (-L0/2., -L0/2.);

  periodic (right);
//  periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-10;
  
  //for (N = 16; N <= 64; N *= 2)
	N=128;
    run();
}

#define WIDTH 0.5
#define EPS 1e-14

event init (t = 0) {
/**
	Initialize microswimmer particles. */
	microswimmers = init_tp_circle(NP);
/**
	If I've got one single particle, set it in the center of the domain.*/
	if(NP==1)
		foreach_particle(){
			p().x = 0.;
			p().y = 0.;
		}
/**
	Initialize particle size. */
	foreach_particle_in(microswimmers){
		p().r = RADIUS;
	}

	compute_SP(microswimmers);
	compute_variable_viscosity();
  mu = muv;

  /**
  The channel geometry is defined using Constructive Solid Geometry. */  
  mask (y > +HEIGHT/2. ? top : none);
  mask (y < -HEIGHT/2. ? bottom : none);
	
	u.n[top]    = dirichlet(0.);
	u.t[top]    = dirichlet(y);
	u.n[bottom] = dirichlet(0.);
	u.t[bottom] = dirichlet(y);
	u.t[right]  = dirichlet(0.); 
	u.n[right]  = dirichlet(0.);
	u.n[left]   = dirichlet(0.);
	u.t[left]   = dirichlet(0.);
}


/**
We check for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
	foreach_particle_in(microswimmers){
		fprintf(stderr,"PARTICLE %+6.5e \n",p().omega);
	}
}

event acceleration (i++){
//	compute_microswimmer_forcing(microswimmers);
//	a = av;
}

event properties (i++)
{
	foreach_particle_in(microswimmers)
		p().r = RADIUS;
	compute_SP(microswimmers);
	compute_variable_viscosity();
  mu = muv;
	vorticity(u,omega);
}
