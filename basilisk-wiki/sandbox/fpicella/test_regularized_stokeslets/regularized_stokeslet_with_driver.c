/*
A clone of the sedimentation_driver.c case
but important because it is a first implementation
capable of multiple-particles.
Here the goal is to add the presence of body forces
in the form of BLOBS
Like [regularized stokeslets](https://epubs.siam.org/doi/10.1137/S106482750038146X)
*/

#include "ghigo/src/myembed.h"
#include "ghigo/src/mycentered.h"
#include "view.h"

#define RADIUS 1.0
#define SEDIMENTATION 0.
#define NPARTICLES 1

#include "fpicella/src/driver-myembed-particles.h"
double THRUST = 100.;
#include "fpicella/src/driver-regularized-stokeslets.h"

/*
Definition of the case and the computational domain
*/
//double HEIGHT = 100.0+1e-2; // container height
//#define container difference(-y+HEIGHT/2.,-y-HEIGHT/2.)
#define container -(sq(x)+sq(y)-sq(L0/2.*0.9))
/*
### Setup
We need a field for viscosity so that the embedded boundary metric
can be accounted. */

face vector muv[];

/*
Define mesh adaptation parameters.*/
int lmin = 4;
int lmax = 6;
#define cmax (1.e-4) // Absolute refinement criteria for the velocity field

FILE * output_file; // global output file pointer

double thetaP = -M_PI/2.;

int main ()
{
	display_control(THRUST,-100,100);
	display_control(thetaP,-10,10);
/*
Space and time are made dimensionless. This is necessary to be able
to use the 'mu = fm' trick.
*/


  size (1. [0]);

  L0 = 16.;

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
		solid(cs,fs,container);//compute the fraction associated to the confinement
		compute_particle_fractions();
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
                        maxlevel = (lmax), minlevel = (lmin));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif

/*
Set the viscosity field in the event *properties*. */
	mu = muv;
/*
Set the acceleration field in the event *acceleration*. */
	a = av;

/*
Boundary condition on non-embedded boundaries */
	u.n[top] = dirichlet(0.);
	u.t[top] = dirichlet(0.);
	u.n[bottom] = dirichlet(0.);
	u.t[bottom] = dirichlet(0.);
	u.n[right] = dirichlet(0.);
	u.t[right] = dirichlet(0.);
	u.n[left] = dirichlet(0.);
	u.t[left] = dirichlet(0.);
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
	foreach_particle()
		fprintf(stderr,"FORCES %+6.5e %+6.5e %+6.5e\n",p().F.x,p().F.y,sqrt(sq(p().F.x)+sq(p().F.y)));
/*
Compute velocities so to guarantee the forece-free condition */
//	velocity_for_force_free();
/*
term 2.*THRUST is places at the center of the embedded body,
it will be opposed to the force applied on Cis-Trans flagella
(represented by the two separate Regularized Stokeslets of 
intesity THRUST) */
	velocity_force_free_sedimentation_propulsion();
///*
//Update particle's location */
//	if(i>10){
//		particle_location_update();
//	}
}

event properties(i++){
	foreach_particle()
		p().theta.z=thetaP;
}

/*
Re-compute and adapt mesh, accounting for embedded boundary as well */
event adapt (i++)
{
	solid(cs,fs,container);//compute the fraction associated to the confinement
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


event logfile (t+=0.1; i<=1000){
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
cylinder of radius 0.25 placed in a container of height 1.0, at center
is around 100 units.
Then, If I set sedimentation force equal to SEDIMENTATION, translational
velocity of a free-sedimenting particle will be around SEDIMENTATION/100.
voilàaaaa
*/
