/**
# Gas bubbles impacting a free-interface and generating a liquid droplet

A 1-mm diameter air bubble rising in a quiescent bath of water,
impacting against the free-surface. For the set of parameters a water droplet is generated and ejected in the gas phase.
Strongly inspired from [bubble example](https://basilisk.fr/src/test/axi_rising_bubble.c) 
(2D-axisymmetric) and the [bubble](https://basilisk.fr/src/examples/bubble.c) one.
*/

/**
![A gas bubble rising towards a liquid interface and generating a bubble.](axi_rising_bubble_to_droplet/movie.mp4)(width="100%")
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"

double bubble_radius = 0.5;
double box_size = 10.;
double endTime = 2.5; 
	
int MAXLEVEL;


int main (int argc, char **argv)
{
  size (box_size);

  MAXLEVEL = 8;
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
			 	- (sq(bubble_radius) - sq(x - box_size*0.2) - sq(y))
				,-x+box_size/2.+1e-4)
				);
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u},(double[]) {0.001,0.001,0.001,0.001},
		 maxlevel = MAXLEVEL);
}
#endif

event movie(i+=1,t<=endTime){
  view(fov=0, tx = 0, ty = -0.5, psi=-pi/2., width = 1200, height = 1200);
	draw_vof ("f", lc = {1, 0, 0}, lw = 2);
	squares ("f", linear = true);
	cells (lc = {0.5, 0.5, 0.5});
	// mirror part of the plot, for ease of view...
	mirror(n={0,1,0}){
	draw_vof ("f", lc = {1, 0, 0}, lw = 2);
	vectors ("u", scale = 0.005);
	};
  save("movie.mp4");
}

// to be validated quantitatively with experiments.
