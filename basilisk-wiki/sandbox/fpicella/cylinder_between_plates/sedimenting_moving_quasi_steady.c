/**
	Single cylindrical particle between parallel plates.
	stemming from [here](http://basilisk.fr/sandbox/fpicella/cylinder_between_plates/sedimentation_driver_moving.c)
*/
double RADIUS = 0.25; // cylinder radius
double HEIGHT = 1.0; // channel height
double nu     = 1.0; // fluid viscosity;

#include "grid/quadtree.h"
#include "ghigo/src/myembed.h"
#include "ghigo/src/mycentered.h"
#include "fpicella/src/driver-myembed-particles.h"
#include "view.h"

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

int lmin = 4;  // Min mesh refinement level (l=6 is 2pt/d)
int lmax = 6; // Max mesh refinement level (l=10 is 36pt/d, l=13 is 292pt/d)
#define cmax (1.e-6) // Absolute refinement criteria for the velocity field

//#define ADAPT 

//#define FLOW (-6.*y*y+1.5) // POISEUILLE
#define FLOW 0.

/**
	### Time-tricker
	This one is the (pseudo-) time step that is used in the Stokes solver.
	It is nothing but a relaxation parameter. I want to set it:
	- low for stability,
	- high for rapid convergence.
	You'll have to work it out to your pleasure. */
#define DT_STOKES 1e-3 // (pseudo-)timestep for the Stokes solver
#define MAXIT_STOKES 250 // maximum number of subiterations for Stokes solver


FILE * output_file; // global file pointer

/**
	For plotting purposes only: timestepper type.
	1 = Explicit Euler
	2 = Adams Bashforth, two steps. */
int timestepper = 0;

int main ()
{
  display_control(RADIUS,0.00,0.5);
  display_control(lmin,1,9);
  display_control(lmax,1,9);
  /**

  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (1. [0]);

  L0 = 2.;

  X0 = Y0 = -L0/2.;

  /**
  We turn off the advection term. The choice of the maximum timestep
  and of the tolerance on the Poisson and viscous solves is not
  trivial. This was adjusted by trial and error to minimize (possibly)
  splitting errors and optimize convergence speed. */

  stokes = true;
  DT = 0.01 [0];//2e-5 [0];

  /**
  We set the tolerance of the Poisson solver. */
  TOLERANCE = 1.e-4;

  /**
  We initialize the grid. */

  N = 1 << (lmax);//initialize at maximum refinement, so to have _child_

  init_grid (N);

  output_file = fopen ("Free_Particle.dat", "w");

	//for(lmax = 7; lmax>=7; lmax -=1){
	////	for(yShift = 0.0; yShift<= 0.2125; yShift += 0.0125)
	////	{
	////		RADIUS = 0.25;
  ////		run();
	////	}

	//	for(RADIUS = 0.1; RADIUS <= 0.4; RADIUS += 0.05)
	//	{
	//		//yShift = 0.;
  //		run();
	//	}
	//}

	for (timestepper = 1; timestepper <=2; timestepper += 1) 
	{
		run();
	}

	fclose(output_file);
}


event init(i=0)
{
  /**
	Initialize particle structures and definitions */
	initialize_particles();
	
	/**
	Setup custom quantities. */
	foreach_particle(){
		p().B.y = -100.;
	}
	/**
	Initialize cs and fs fields associated to the particles */
	compute_particle_fractions();

//	/**
//	Horizontally-oriented channel */
//  mask (y > +HEIGHT/2. ? top : none);
//  mask (y < -HEIGHT/2. ? bottom : none);
	/**
	Vertically-oriented channel */
  mask (x > +HEIGHT/2. ? right : none);
  mask (x < -HEIGHT/2. ? left : none);

#ifdef ADAPT
#if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */

  astats ss;
  int ic = 0;
  do {
    ic++;
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
                        maxlevel = (lmax), minlevel = (lmin));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
#endif
  /**
  We set the viscosity field in the event *properties*. */

  mu = muv;
	
	/**
	Set boundary conditions on non-embed bodies */

  u.n[top]    = dirichlet(0.);
  u.t[top]    = dirichlet(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);
  u.t[right]  = dirichlet(0.);
  u.n[right]  = dirichlet(FLOW);
  u.n[left]   = dirichlet(FLOW);
  u.t[left]   = dirichlet(0.);
#if dimension == 3
	u.n[front]  = dirichlet(0.);
	u.t[front]  = dirichlet(0.);
	u.n[back]   = dirichlet(0.);
	u.t[back]   = dirichlet(0.);
#endif

  
/**
	Boundary condition on the particles (or any other embedded boundary... ) */
	u.n[embed]  = dirichlet(velocity_noslip_x(point,x,y,z));
	u.t[embed]  = dirichlet(velocity_noslip_y(point,x,y,z));
	uf.n[embed] = dirichlet(velocity_noslip_x(point,x,y,z));
	uf.t[embed] = dirichlet(velocity_noslip_y(point,x,y,z)); // pOmg.x = pOmg.y
}



event properties (i++) // refresh particle's BC on embed at each iteration...!
{
/**
###	Update particle's location.
	I let it move only after few iterations,
	to avoid numerical oscillations and startup timesteppers.  */
	//if(i>10){
	//	particle_location_update();
	//}
/**
	In case I want to change some of the particle's properties in runtime. */
	foreach_particle()
		p().r = RADIUS;
	compute_particle_fractions();

  foreach_face()
    muv.x[] = (nu)*fs.x[];
  boundary ((scalar *) {muv});
/**
	Set velocities inside the particles to zero.
	For cosmetic purposes only. */
	foreach()
		if(cs[]==0.)
			foreach_dimension()
				u.x[] = 0.;
}

#ifdef ADAPT
event adapt (i++) // and not classic adapt, so to keep mesh identical during subiterations...
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (lmin));
}
#endif

/**
	### Stop simulation when particle is impacting the ground */
event logfile (t += dt; i <= 5000) {
	// Dummy output
	foreach_particle(){
		fprintf(stderr,"PARTICLE %+6.5e %+6.5e %+6.5e %+6.5e \n",p().u.x,p().u.y,p().F.x,p().F.y);
		fprintf(output_file,"%04d %04d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e \n",timestepper,N,dt,p().u.x,p().u.y,p().F.x,p().F.y,p().x,p().y);
		fflush(output_file);
		if(p().y<-L0/2.+p().r)
			return 1.;
	}
}


  /**
	# Quasi-steady solution.
	I take advantage of the end_timestep hook, from
	centered.h formulation, so to get quasi-steady solutions.
	I can use it so to converge as well on the computation
	of the particle's velocities required to get the 
	force-torque free condition.
	*/
	/**
	# Time-tricker
	before getting into the solver, I set the _quasi-steady_ timestep to
	an arbitrary low value...in any case, it is a fictitious timestep,
	required only to converge the prediction-correction iterator! */
event begin_timestep(i++){
	dt = DT_STOKES;
}
	
event end_timestep(i++){
	/**
	Re-run prediction-projection steps to converge
	in a single timestep.*/
	scalar un[];
	int j = 0;
	double du = 1.;
	/**
	Store previous velocities...required for
	Adams Bashforth, linear multistep technique.*/
	foreach_particle()
		foreach_dimension()
			p().u2.x = p().u.x;

	while (j<MAXIT_STOKES){
		/**
		Real, hard-core computation:
		solving for the NS problem.*/
		event ("viscous_term");
		event ("acceleration");
		event ("projection");	
		/**
		Lightweight: compute forces on particles.*/
		hydro_forces_torques();
		/**
		Compute the velocity that I would need so to guarantee
		the force free condition.*/
		velocity_for_force_free();

		/**
		Stopping criteria.*/
  	du = change (u.y, un);
		foreach()
			un[] = u.y[];
		j += 1;
//		fprintf(stderr,"ITERATION %03d du %+6.5e\n",j,du);
		/**
			Convergence criterion */
		if(du<1e-4)
			j = MAXIT_STOKES;
	}
	/**
	Time-tricker, reset dt.*/
	dt = DT;

	/**
	Update particle location once that convergence is obtained. */
	if(timestepper == 1)
	particle_location_update();
	if(timestepper == 2){
	/** 
	Initialize with Explicit Euler.*/
	if (i <= 1)
		foreach_particle()
			foreach_dimension()
				p().x += p().u.x*dt;
	else
	/**
	Adams Bashforth, linear multistep. */
		foreach_particle()
			foreach_dimension()
				p().x += p().u.x*dt*3./2. - p().u2.x*dt*1./2.;
	}
}


/**
# Results
~~~gnuplot Sedimenting cylindrical particle
set grid
set xlabel "Particle-wall gap"
set ylabel "Sedimenting velocity"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot 'Free_Particle.dat' using ($1==1 ?$9+0.75 : 1/0):($1==1 ? $5 : 1/0) with points pt 7 ps 2 lc rgb "green" title "Explicit Euler",\
     'Free_Particle.dat' using ($1==2 ?$9+0.75 : 1/0):($1==2 ? $5 : 1/0) with points pt 5 ps 1 lc rgb "red" title "Adams Bashforth"
~~~
For this configuration, I do not observe any significant advantage of AB.
*/
