/**
# Gas bubbles interacting with particles

A 1-mm diameter air bubble rising in a quiescent bath of water.
Strongly inspired from [bubble example](https://basilisk.fr/src/test/axi_rising_bubble.c) 
(2D-axisymmetric) and the [bubble](https://basilisk.fr/src/examples/bubble.c) one.

### Particles
Applying [Antoon's code](https://basilisk.fr/sandbox/Antoonvh/pc.c), I employ the computed velocity field to get the dynamics of non-inertial (Stokes, black, left) and inertial (red, right) particles seeded within the fluid.
I should get an idea of the particle-wake interaction around the spherical bubble.
*/

/**
![A gas bubble rising in a quiescent fluid. Tracer (non-inertial) particles on the left (black), inertial particles on the right (red).](axi_rising_bubble_with_particles/movie.mp4)(width="100%")
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"

/**
I'd like to somehow track tracer and intertial particles within the flow.
To do so, I try to use Antoon's particle implementation.
The inspiration comes from his case of the [2D turbulence](https://basilisk.fr/sandbox/Antoonvh/pc.c).*/
#include "../../Antoonvh/tracer-particles.h"
#include "../../Antoonvh/stokes-particles.h"
#include "../../Antoonvh/scatter2.h"
Particles flow, heavy;
 

double bubble_radius = 0.5;
double box_size = 10.;
double endTime = 5.0; 
	
int MAXLEVEL;


int main (int argc, char **argv)
{
  size (box_size);

  MAXLEVEL = 7;
#if !TREE
  N = 1 << MAXLEVEL;
#endif

/** The case is made dimensionless using the same approach as the one described in the [bubble example](https://basilisk.fr/src/examples/bubble.c). */

#define RHOR 1000.
#define MUR 100.

	// I consider gas bubbles that are so small that they stay almost spherical when rising.
	// according to Saffman 1956, it should have size of almost 1 mm
	// Ga \approx 1000, Bo \approx 0.1.
	// To be checked again, but it seems to work properly

	const double Ga = 1000.; // gravity vs viscous forces
	const double Bo = 0.1;   // gravity vs capillary forces

	rho1 = 1. [0];
  rho2 = rho1/RHOR;
  mu1 = 1./Ga;
  mu2 = 1./(MUR*Ga);
  f.sigma = 1./Bo;

  // density
	G.x = -1.; // acceleration of gravity is set to unity... see examples/bubble.c

	run();
}


event init (t = 0)
{
// Mask upper part of the domain
// that corresponds to the zone that is radially further
// from the system's axis of symmetry.
	mask (y > L0/2 ? top : none);

#if TREE
	// refine bubble location
  refine (sq(2.*bubble_radius) - sq(x - box_size*0.25) - sq(y) > 0 &&
	  level < MAXLEVEL);
	// refine interface location
  refine (-x+box_size/2. && level < MAXLEVEL);
#endif

  fraction (f,intersection(
			 	- (sq(bubble_radius) - sq(x - box_size*0.1) - sq(y))
				//,-x+box_size/2.+1e-4 // free-interface
				,1. // nothing
				)
				);
/*
Initialize Particles.*/
  new_tracer_particles (0);
  new_inertial_particles (0);
/**
Particles are set at every cell center*/
	flow  = init_tp_cells();
  heavy = init_ip_cells();
/**
Now I remove particles that sit within the gas phase.*/
	remove_particles(flow,  interpolate(f, x, y, z) < 0.999);
	remove_particles(heavy, interpolate(f, x, y, z) < 0.999);

/**
And I assigh how _heavy_ my particles are.*/
  foreach_particle_in(heavy)
    p().u2.z = 0.1;//0.001;  //Set relaxation timescale
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u},(double[]) {0.001,0.001,0.001,0.001},
		 maxlevel = MAXLEVEL);
}
#endif

event movie(t+=endTime/100.,t<=endTime){
  view(fov=0, tx = 0, ty = -0.5, psi=-pi/2., width = 1200, height = 1200);
	draw_vof ("f", lc = {0, 0, 1}, lw = 2);
	//squares ("f", linear = true);
	cells (lc = {0.75, 0.75, 0.75});
	scatter (flow , s = 3);
	// mirror part of the plot, for ease of view...
	mirror(n={0,1,0}){
	draw_vof ("f", lc = {0, 0, 1}, lw = 2);
	//vectors ("u", scale = 0.005);
	scatter (heavy, s = 3, pc = {1, 0, 0});
	};
  save("movie.mp4");
}

// to be validated quantitatively with experiments.