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

/**
	Output file, file pointer. */
FILE * output_file; // global file pointer

/**
	Number of particles in simulation. */
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

  output_file = fopen ("Free_Particle.dat", "w");
  for (N = 16; N <= 256; N *= 2)
		for (RADIUS = 0.1; RADIUS <= 0.4; RADIUS += 0.05)
	//N=128;
    run();
	fclose(output_file);
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
scalar un[]; // field that will contain previous solution...
event logfile (t += 0.01; i <= 1000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-6){ // since I'm looking for a steady solution, this must be quite small...
		foreach_particle_in(microswimmers)
  	fprintf(output_file,"%+3.2e %03d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
                L0, N, ETA, p().r, p().x, p().y, p().theta, p().u.x, p().u.y, p().omega);
  	fflush(output_file);
    return 1; /* stop */
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


/**
# Force-Torque free 2D cylinder in plane couette flow
~~~gnuplot Free cylinder in plane couette, Stokes flow
set grid
set xlabel "Radius"
set ylabel "Omega"
set title "Free cylinder in plane Couette"
set key right left

plot \
    'Free_Particle.dat' using ( $2==256 ? $4 : 1/0 ):( $2==256 ? +$10 : 1/0 ) with points pt 13 ps 2 lc rgb "#fb9f3a" title "N = 256",\
    'Free_Particle.dat' using ( $2==128 ? $4 : 1/0 ):( $2==128 ? +$10 : 1/0 ) with points pt 9 ps 2 lc rgb "#b63679" title "N = 128",\
    'Free_Particle.dat' using ( $2==64 ? $4 : 1/0 ):( $2==64 ? +$10 : 1/0 ) with points pt 5 ps 2 lc rgb "#3b0f70" title "N = 64",\
    'Free_Particle.dat' using ( $2==32 ? $4 : 1/0 ):( $2==32 ? +$10 : 1/0 ) with points pt 7 ps 2 lc rgb "#000004" title "N = 32",\
    '../../cylinder_between_plates/Dvinsky_Popel_1987_fig_11.dat' using ( $1):( $2) \
     with lines lw 5 lc rgb "black" title "Dvinsky Popel 1987, fig 11"

~~~
*/
