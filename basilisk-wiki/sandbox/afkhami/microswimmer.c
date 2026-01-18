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

  L0 = 2.;

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

  //TOLERANCE    = 1.e-6;
  //TOLERANCE_MU = 1.e-6;

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

	run();

	fclose(output_file);
}


event init(i=0)
{
  /**
	Initialize particle structures and definitions */
	initialize_particles();
	foreach_particle(){
		p().B.y = -sq(p().r)*100.;
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
###	Update particle's location.
	I let it move only after few iterations,
	to avoid numerical oscillations and startup timesteppers.  */
	if(i>10){
		particle_location_update();
	}
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
	### Stop simulation when particle is impacting the ground */
event logfile (t += 0.01; i <= 5000) {
	// Dummy output
	foreach_particle()
		fprintf(stderr,"PARTICLE %+6.5e %+6.5e %+6.5e %+6.5e \n",p().u.x,p().u.y,p().F.x,p().F.y);
	foreach_particle()
		if(p().y<=-L0/2.+p().r){
			return 1.;
		}
}

/**
# Results
~~~gnuplot Sedimenting cylindrical particle
set grid
set xlabel "particle-wall gap"
set ylabel "Normalized sedimenting velocity"
set yrange [0:1.5]


# Looking at Dvinsky Popel 1986, fig 7a, I can get what would be the sedimenting
# velocity of my particle. For a particle of radius 0.25, I should get
# around U_sedimentation = 6.300 e-4.
# Now, I use a sedimentation force that is set here (line 100 in this file)
# to be 100 times the one of Dvinsky Popel. But since I'm in a linear regime
# I should just get a velocity that is 100 times bigger.
# then, my target velocity is for a cylindrical particle settling between 2 plates
# but still FAR ENOUGH from the top-bottom, is U_far=6.3e-2.
# I will use this value to normalize my observed velocity.

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot 'particle_000.dat' using ($3+1-0.25):(-$6/6.3e-2) \
     with points pt 7 ps 2 lc rgb "green" title "R=0.25"
~~~

### Influence of hydrodynamic forces on sedimenting velocity
- Velocity of a confined cylinder, but far from container's end
	it matches the prediction of Dvinsky Popel 1986
- Velocity computed as a force-free condition on embedded boundaries
- Location is advanced using a simple Explicit Euler timestepper.
- Velocity is normalized wrt the one predicted by Dvinsky Popel 1987
	figure 7a.
- Then, for large gap, I should get strictly 1...and it is!
### Particle slows-down as it is approaching the container's wall
- Full hydrodynamic interaction. 
- No additional potential is added (Like Lennard - Jones...)
- I can not capture lubrication ( I should get infinite force->zero velocity
	when the gap goes to zero), but still, there's the qualitative behaviour.
*/

/*
Display a solution...*/

event snapshot (i += 5)
{
  view (fov = 20, camera = "front",
	tx = 0., ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  squares ("u.y");
	vectors ("u", scale = 0.1);
  save ("output.mp4");
}

/**
# Moving, sedimenting cylinder
![Sedimentation of a 2D confined cylinder, getting closer to a wall. Stokes regime](sedimentation_driver_moving/output.mp4)(loop)
*/
