/**
# Single tHree micro-swimmer

tHree
*/

#include "navier-stokes/centered.h"
#include "view.h"
/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */

double HEIGHT = 32.0;
double THETA = M_PI/2.*1;
double RADIUS = 1.0;
double DISTANCE = 4.0;
double FACTOR = 1.;
double B = 1.; // reorientation owing to gravity
double BUOYANCY = 0.;
double THRUST   = 4.;
/** Propulsive system's geometry.*/
/** monodisperse family, being all identical (all puller, pusher, or whatever...*/
double beta = 1.0;
double dOverR = 2.0;

#include "fpicella/src/periodic-shift-treatment.h"

#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double theta; double thetan; double omega; double omega2; double r; \
coord B; double Thrust;

//#define alpha +2.0
//#define beta  +2.0

//#define Particle_Advection_ON 1
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

  size (256 [0]);
  //DT = HUGE [0];
  DT = 100.0 [0];// higher DT, equivalent to higher viscosity
									// it's kind of the stiffness of the problem
									// but since I'm looking for a steady solution
									// I can cranck it to high level
									// so to obtain a solution more rapidly.
  
  origin (-L0/2., -L0/2.);

  //periodic (right);
  //periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-3;

  output_file = fopen ("Output.dat", "w");

  
//  for (N = 32; N <= 256; N *= 2)
//		for(RADIUS = 0.05; RADIUS <= 0.45; RADIUS += 0.05)
	N=256;
	for(dOverR = L0/2.; dOverR >= RADIUS*0.1; dOverR *= 0.9)
  	run();

	fclose(output_file);
}

#define WIDTH 0.5
#define EPS 1e-14

event init (t = 0) {
/**
	Initialize microswimmer particles. */
	//microswimmers = init_tp_circle(NP);
	microswimmers = init_tp_square(NP);
	foreach_particle_in(microswimmers){
		p().B.y=BUOYANCY;
		p().r = RADIUS;
		p().Thrust = THRUST;
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


event acceleration (i++){
	compute_microswimmer_forcing_smooth(microswimmers);
//	compute_repulsion(microswimmers,3.0,false,false,false,false,12.); // look at the driver-tHree file
																																		// to know what argument is meant
																																		// to do ;)
	a = av; // assign acceleration term. That's capital!
}

event properties(i++){
	vorticity(u,omega);
	compute_bodyPlot(microswimmers);
}

event logfile ( i += 1; i <= 2500) {
  double du = change (u.y, un);
  if (i > 5 && du < 1e-5){ // since I'm looking for a steady solution, this must be quite small...
		foreach_particle()
	  fprintf(output_file,"%+3.2e %+04i %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e \n",
            L0, N, RADIUS, beta, dOverR, p().Thrust, BUOYANCY, p().u.y);
		fflush(output_file);

    return 1; /* stop */
  }
}

/**
### Proper definition of _stopping_ event for the solver.*/

/**
~~~pythonplot tHree swimming velocity, as a function of the geometry of the propulsion system.
import numpy as np
import matplotlib.pyplot as plt

cmap = plt.get_cmap('plasma')

data = np.loadtxt("Output.dat")

fig, ax = plt.subplots(figsize=(4,3))

ax.plot(data[:,4], data[:, 7]/data[0,7], 'o', color=cmap(0.25), markersize=8, label="2D, confined, tHree")
dR= np.linspace(1.0,128,100); aOd = 1/dR;
alpha = np.arctan(1.0);
V0=(1+(-3/2*aOd+1/2*aOd**3)*(np.sin(alpha))**2+(-3/4*aOd-1/4*aOd**3)*(np.cos(alpha)**2));
ax.plot(dR,V0, color=cmap(0.0),lw=5, label="3D, unbounded, analytical");

### Scaling laws...
L1 = aOd;
ax.plot(dR,1-aOd,'--', color=cmap(1.0), label = "$ 1-R/d$");
b=1;
#L2 = np.log(1+b*dR)/np.log(1+b*dR[-1]);
#ax.plot(dR,L2,':k', label= "$\propto ln(1+d/R)/ln(1+d_{max}/R)$");
L2 = 1-np.exp(-dR*0.065);
ax.plot(dR,L2,'-', color = cmap(0.75), label= "$1-exp(-d/R*k)$");
plt.xlabel('$d/R$');
plt.ylabel('$U_{swim}$');
ax.legend()
plt.tight_layout();
plt.savefig('single_tHree_velocity_as_geometry_PALOTAI.svg')
~~~
*/
