/**
	Single cylindrical particle between parallel plates.
*/
double RADIUS = 0.25; // cylinder radius
double HEIGHT = 1.0; // channel height
double nu     = 1.0; // fluid viscosity;

//double SEDIMENTATION = RADIUS*RADIUS; //
//#define SEDIMENTATION RADIUS*RADIUS

#include "grid/quadtree.h"
#include "ghigo/src/myembed.h"
#include "ghigo/src/mycentered.h"
#include "fpicella/src/driver-myembed-particles.h"

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

int lmin = 4;  // Min mesh refinement level (l=6 is 2pt/d)
int lmax = 7; // Max mesh refinement level (l=10 is 36pt/d, l=13 is 292pt/d)
#define cmax (1.e-6) // Absolute refinement criteria for the velocity field

FILE * output_file; // global file pointer

int main ()
{
  display_control(RADIUS,0.00,0.5);
  display_control(lmin,1,9);
  display_control(lmax,1,9);
	//display_control(SEDIMENTATION,-1000,1000);
  /**

  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (1. [0]);

  L0 = 4.;

  X0 = Y0 = -L0/2.;

  /**
  We turn off the advection term. The choice of the maximum timestep
  and of the tolerance on the Poisson and viscous solves is not
  trivial. This was adjusted by trial and error to minimize (possibly)
  splitting errors and optimize convergence speed. */

  stokes = true;
  DT = 1e0;//2e-5 [0];

  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6;

  /**
  We initialize the grid. */

  N = 1 << (lmax);//initialize at maximum refinement, so to have _child_

  init_grid (N);

  output_file = fopen ("Free_Particle.dat", "w");

	for(lmax = 7; lmax>=7; lmax -=1){
	//	for(yShift = 0.0; yShift<= 0.2125; yShift += 0.0125)
	//	{
	//		RADIUS = 0.25;
  //		run();
	//	}

		for(RADIUS = 0.1; RADIUS <= 0.4; RADIUS += 0.05)
		{
			//yShift = 0.;
  		run();
		}
	}

//	run();

	fclose(output_file);
}


event init(i=0)
{
  /**
	Initialize particle structures and definitions */
	initialize_particles();
	foreach_particle(){
		p().B.y = -sq(p().r);
//		p().omega.x = +0.001;
//		p().omega.y = +0.001;
//		p().omega.z = +0.001;
	}
	/**
	Initialize cs and fs fields associated to the particles */
	compute_particle_fractions();
	/**
	Employ masks to identify the presence of container's wall.*/
//	/**
//	Horizontally-oriented channel */
//  mask (y > +HEIGHT/2. ? top : none);
//  mask (y < -HEIGHT/2. ? bottom : none);
	/**
	Vertically-oriented channel */
  mask (x > +HEIGHT/2. ? right : none);
  mask (x < -HEIGHT/2. ? left : none);
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
  /**
  We set the viscosity field in the event *properties*. */

  mu = muv;
	
	/**
	Set boundary conditions */
  
  u.n[left]  = dirichlet(0.);
  u.t[left]  = dirichlet(0.);
    p[left]  = neumann(0.);
  u.n[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
    p[right] = neumann(0.);

  u.n[top]  = dirichlet(0.);
  u.t[top]  = dirichlet(0.);
    p[top]  = neumann(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);
    p[bottom] = neumann(0.);
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
	In case I want to change some of the particle's properties in runtime. */
	foreach_particle()
		p().r = RADIUS;
	compute_particle_fractions();
/**
### Compute particle's dynamics.*/
	hydro_forces_torques();
	velocity_for_force_free();

  foreach_face()
    muv.x[] = (nu)*fs.x[];
  boundary ((scalar *) {muv});
}

/**
### Compute particle's dynamics. Exploit the _acceleration_ event*/
event acceleration (i++){
	


}

event adapt (i++) // and not classic adapt, so to keep mesh identical during subiterations...
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (lmin));
/**
	Set velocities inside the particles to zero.
	For cosmetic purposes only. */
	foreach()
		if(cs[]==0.)
			foreach_dimension()
				u.x[] = 0.;
}

/**
We look for a stationary solution. */
scalar un[]; // field that will contain previous solution...
event logfile (t += 0.01; i <= 1000) {
	// Dummy output
	foreach_particle()
		fprintf(stderr,"PARTICLE %+6.5e %+6.5e %+6.5e %+6.5e \n",p().u.x,p().u.y,p().F.x,p().F.y);
  double du = change (u.y, un);
  if (i > 0 && du < 1e-7){ // since I'm looking for a steady solution, this must be quite small...
//                           // making the simulation quite expensive
		foreach_particle()
  	fprintf(output_file,"%+3.2e %02d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
                L0, lmax, RADIUS,  0., p().u.x, p().u.y, p().omega.x,p().F.x,p().F.y,p().T.x);
  	fflush(output_file);
    return 1; /* stop */
  }
}

/**
# Forces on the cylinder
~~~gnuplot Sedimenting cylindrical particle
set grid
set xlabel "Radius"
set ylabel "Sedimenting Velocity"
set title "Sedimenting cylinder, y direction, embed + mask + driver"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot '../Dvinsky_Popel_1986_fig_07a.dat' using ( $1):( $2) \
     with lines lw 5 lc rgb "red" title "Dvinsky Popel 1986, fig 07 a",\
     'Free_Particle.dat' using ( $4==0 && $2==7 ? $3 : 1/0 ):( $4==0  && $2==7 ? -$6 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "blue" title "lmax 7",\
     'Free_Particle.dat' using ( $4==0 && $2==8 ? $3 : 1/0 ):( $4==0  && $2==8 ? -$6 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "green" title "lmax 8"
~~~
Dead-on spot, a nice result!
*/


