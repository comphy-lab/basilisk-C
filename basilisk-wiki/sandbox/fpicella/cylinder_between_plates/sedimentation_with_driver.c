/*
A clone of the sedimentation.c case
but important because it is a first implementation
capable of multiple-particles.
*/

#include "ghigo/src/myembed.h"
#include "ghigo/src/mycentered.h"
#include "view.h"

#define RADIUS 0.25
#define SEDIMENTATION 25.
#define NPARTICLES 1

#include "fpicella/src/driver-myembed-particles.h"

/*
Definition of the case and the computational domain
*/
double HEIGHT = 1.0+1e-2; // channel height
#define channel difference(-y+HEIGHT/2.,-y-HEIGHT/2.)
/*
### Setup
We need a field for viscosity so that the embedded boundary metric
can be accounted. */

face vector muv[];
/*
Define mesh adaptation parameters.*/
int lmin = 4;
int lmax = 7;
#define cmax (1.e-6) // Absolute refinement criteria for the velocity field

FILE * output_file; // global output file pointer

int main ()
{
/*
Space and time are made dimensionless. This is necessary to be able
to use the 'mu = fm' trick.
*/


  size (1. [0]);

  L0 = 4.;

  X0 = Y0 = -L0/2.;

/*
Turn off the adcection term, set DT and tolerances for Poisson solver.*/
  stokes = true;
  DT = 1e0 [0];//2e-5 [0];
	TOLERANCE = 1.e-3; // with extra-low tolerance, it runs also if pressure is not well computed
											// I should adress this point, sometime...
/* Initialize grid */
	N = 1 << (lmax);
	init_grid(N);

	run();
}

/* Initialize */
event init (i=0)
{
/*Initialize particle list, structure...*/
	initialize_particles();

/*When using TREE, refine mesh around the embedded boundary*/
#if TREE
  astats ss;
  int ic = 0;
  do {
    ic++;
		solid(cs,fs,channel);//compute the fraction associated to the confinement
		compute_particle_fractions();
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
                        maxlevel = (lmax), minlevel = (lmin));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif

/*
Set the viscosity field in the event *properties*. */
	mu = muv;

/* Set up sedimentation forces as well... */
	foreach_particle()
		p().B.x = SEDIMENTATION;

}



event properties(i++){
	foreach_face()
		muv.x[] = fs.x[];
	boundary ((scalar *) {muv});
/*
Boundary conditions on embedded boundaries.
This is done automatically with respect to Rigid Body Motion
of each particle, as well as for the domain's container. */
	u.n[embed]  =  dirichlet(velocity_noslip_x(point,x,y,z));
	u.t[embed]  =  dirichlet(velocity_noslip_y(point,x,y,z));
	uf.n[embed] = dirichlet(velocity_noslip_x(point,x,y,z));
	uf.n[embed] = dirichlet(velocity_noslip_y(point,x,y,z));
	p[embed]  = neumann(0.);
	pf[embed] = neumann(0.);
/*
Compute forces on particles */
	hydro_forces_torques();
/*
Compute velocities so to guarantee the forece-free condition */
	velocity_for_force_free();
/*
Update particle's location */
	if(i>10){
		particle_location_update();
	}
}

/*
Re-compute and adapt mesh, accounting for embedded boundary as well */
event adapt (i++)
{
	solid(cs,fs,channel);//compute the fraction associated to the confinement
	compute_particle_fractions();
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (lmin));
	/*
	Cleanup velocity
	for visualisation purposes only
	*/
	foreach()
		foreach_dimension()
			u.x[]*=cs[];
}


event logfile (t+=0.01; i<=100){
}

event snapshot (i += 1)
{
  view (fov = 20, camera = "front",
	tx = 0., ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  squares ("u.x", min = 0);
  save ("mesh.mp4");
}

/**
# Moving, sedimenting cylinder
![Mesh](sedimentation_with_driver/mesh.mp4)(loop)
~~~gnuplot Sedimenting cylindrical particle
set grid
set xlabel "Position"
set ylabel "Velocity"
set title "Sedimenting cylinder"

set logscale y

plot \
    'particle_000.dat' using ($2):($5) with points pt 7 ps 2 lc rgb "green" title "computed", \
    '-' with lines lt 1 lw 4 lc rgb "red" title "Prediction from Dvinsky Popel 1986"
0.0 0.25
0.2 0.25
e
~~~
Why this prediction?
From plot 2a in Dvinsky Popel 1986, we know that the force on a 
cylinder of radius 0.25 placed in a channel of height 1.0, at center
is around 100 units.
Then, If I set sedimentation force equal to SEDIMENTATION, translational
velocity of a free-sedimenting particle will be around SEDIMENTATION/100.
voilàaaaa
*/
