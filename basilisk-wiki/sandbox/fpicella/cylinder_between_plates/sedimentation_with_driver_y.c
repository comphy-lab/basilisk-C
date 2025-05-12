/**
# Sedimenting passive cylinder in a channel
	2D, Steady and _Stokes_ flow.
	Example taken from (Dvinsky Popel 1986)[https://doi.org/10.1016/0045-7930(87)90031-4]
*/

//#include "embed.h"
//#include "navier-stokes/centered.h"

#include "ghigo/src/myembed.h"
#include "ghigo/src/mycentered.h"
#include "view.h"

#define RADIUS 0.25
#define SEDIMENTATION 1000.0
#define NPARTICLES 1
/*
Turned off, so that it is not taken into account within velocity_noslip*/
//#define imposedFlow 1
coord imposedU = {0.0,0.0,0.0}; // externally imposed flow
#include "fpicella/src/driver-myembed-particles.h"

/**
A flag for mesh adaptation... */
#define MeshAdaptation 1

/*
Iterative solving...
*/
double radius = 0.25;
double xShift = 0.;
double yShift = 0.;


/*
Mesh adaptation parameters */
int lmin = 3;  // Min mesh refinement level (l=6 is 2pt/d)
int lmax = 6; // Max mesh refinement level (l=10 is 36pt/d, l=13 is 292pt/d)
#define cmax (1.e-4) // Absolute refinement criteria for the velocity field

FILE * output_file; // global file pointer

int main()
{
  output_file = fopen ("Particle_Forces.dat", "w");
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (4. [0]);
  //DT = HUGE [0];
	/**
	DT = HUGE will not work with my simple timesteppers! */
  DT = 1.0 [0];
  
  origin (-L0/2., -L0/2.);
//  periodic (right);
//  periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-3;
  
  //for (N = 16; N <= 64; N *= 2)
  N = 1 << (lmax);//initialize at maximum refinement, so to have _child_
  init_grid (N);
	//N = 128;

//    run();


//	for(lmax = 9; lmax>=7; lmax -=1){
//		for(yShift = 0.0; yShift<= 0.2125; yShift += 0.0125)
//		{
//			RADIUS = 0.25;
//  		run();
//		}

//		for(radius = 0.1; radius <= 0.4; radius += 0.05)
//		{
////			yShift = 0.;
//  		run();
//		}
////	}

	run();

	fclose(output_file);
}

scalar un[];

#define HEIGHT 1.0
#define EPS 1e-14
#define channel difference(-x+HEIGHT/2.+EPS,-x-HEIGHT/2.-EPS) // y-aligned channel
//#define channel difference(-y+HEIGHT/2.+EPS,-y-HEIGHT/2.-EPS) // x-aligned channel

event init (t = 0) {
/**
Initialize particle list, structure...*/
	initialize_particles();
/**
	Set particle properties */
	foreach_particle()
		p().B.y = -SEDIMENTATION*sq(p().r);

  /**
  Viscosity is unity. */
  
  mu = fm;

  /**
  The channel geometry is defined using Constructive Solid Geometry. */  
	solid(cs,fs,channel);
 /**
	Add as well the presence of the particle. */
	compute_particle_fractions();

#ifdef MeshAdaptation
/* 
	Automatic mesh adaptation */
#if TREE
  astats ss;
  int ic = 0;
  do {
    ic++;
		solid(cs,fs,channel);
		compute_particle_fractions();
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
                        maxlevel = (lmax), minlevel = (lmin));
  } while ((ss.nf || ss.nc) && ic < 100);

#endif
#endif

  /**
  The boundary condition is zero velocity on the embedded boundaries. */

  u.n[embed] = dirichlet(velocity_noslip_x(point,x,y,z));
  u.t[embed] = dirichlet(velocity_noslip_y(point,x,y,z));
/*
	myembed.h requires to define uf fields as well... */
  uf.n[embed] = dirichlet(velocity_noslip_x(point,x,y,z));
  uf.t[embed] = dirichlet(velocity_noslip_y(point,x,y,z));

  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */

  for (scalar s in {u})
    s.third = true;
  
}

event particle_evolution(i++){
/**
	### Compute forces on particles. */
	hydro_forces_torques();
/**
	### Determine velocity so to be force free.*/
	velocity_for_force_free();
/**
	### Advect particle location 
	A simple way to avoid instabilities, is to let the solver
	adapt by itself for a couple of iterations...
	*/
	if(i>5){ 
		particle_location_update();
	}
	
}

#if MeshAdaptation
event adapt (i++) // and not classic adapt, so to keep mesh identical during subiterations...
{
		solid(cs,fs,channel);
		compute_particle_fractions();
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (lmin));
	/**
	Cleanup fractions */
	foreach()
		foreach_dimension()
			u.x[] *=cs[];
}
#endif

/**
We check for a stationary solution. */

event logfile (t += 0.01; i <= 100) {
}

/*
Display a solution...*/

event snapshot (i += 1)
{
  view (fov = 20, camera = "front",
	tx = 0., ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  squares ("u.y", min = 0);
  save ("output.mp4");
}

/**
# Moving, sedimenting cylinder
![Mesh](sedimentation_with_driver_y/output.mp4)(loop)

#Need a quantitative plot, but it's ok...
*/
