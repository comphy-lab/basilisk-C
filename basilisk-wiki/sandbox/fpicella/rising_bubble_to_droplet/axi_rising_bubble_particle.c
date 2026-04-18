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
#include "tracer-particles.h"
#include "stokes-particles.h"
#include "scatter2.h"
Particles flow, heavy;
 
/**
Domain size, final simulation time, timestepping for post-processing purposes.*/
// BUBBLE AVERAGE DIAMETER IS SET ALWAYS TO 1!
double box_size = 40.;
double endTime = 20.0; 
double stepTime= 0.25;

/**
The only physical parameter I can change in my problem. Average bubble size (i.e. bubble volume...)
injected within the system.*/
double dBubble; // global variable, what is the present bubble diameter for which I'm doing my test?
/**
Phisical parameters of the problem. Air and water.*/
double rhoOut = 1000.0;
double nuOut  = 1.0e-6;
double sigma  = 0.072;
double gCons  = 9.81;
	
int MAXLEVEL;

/**
# First-order analytical bubble shape (Bond number controlled)

We approximate the static bubble shape as a first-order perturbation of a
sphere in the Bond number

$$
Bo = \frac{\Delta \rho g R^2}{\gamma}.
$$

- `Bo = 0` → exact sphere  
- valid for moderate Bond numbers  
- does not capture necking or detachment
*/

static inline double P2 (double mu)
{
  return 0.5*(3.*sq(mu) - 1.);
}

static inline double bubble_phi_axi_bary_unit (double x, double y,
                                               double xBary,
                                               double Bo)
{
  const double eps = Bo/6.;
  const double vf = 1. + (3./5.)*sq(eps) + (2./35.)*eps*sq(eps);

  // radius corrected so that the bubble volume equals that of a sphere of diameter 1
  const double Rs = 0.5*pow(vf, -1./3.);

  const double X = x - xBary;
  const double R = sqrt(sq(X) + sq(y));

  if (R < 1e-14)
    return -(Rs*(1. + eps));

  return -(Rs*(1. + eps*P2(X/R)) - R);
}


int main (int argc, char **argv)
{
  size (box_size);

	origin(-5.,0.); // domain shifted so not to interfere with bubble.

  MAXLEVEL = 8;
#if !TREE
  N = 1 << MAXLEVEL-2;
#endif

/** The case is made dimensionless using the same approach as the one described in the [bubble example](https://basilisk.fr/src/examples/bubble.c). */

#define RHOR 1000.
#define MUR 100.

/**
For my convenience, I define Ga and Bo as macros... */
#define Ga (sqrt(gCons*(dBubble)*(dBubble)*(dBubble)*1e-9)/(nuOut))
#define Bo (rhoOut*(gCons)*(dBubble)*(dBubble)*1e-6/(sigma))

/**
I the only variable parameter here is the bubble diameter dBubb.
I provide it naively within a list .
As a convention, here dBubble is provided in millimiter.
This is the only _dimensional_ parameter, on which I will work out all the other
dimensionless parameters (Bond and Galillei) that drives the problem.
*/
	//double dBubble_list[] = {0.5,1.0,1.5,2.0,2.5};
	double dBubble_list[] = {1.0};
  int n = sizeof(dBubble_list)/sizeof(dBubble_list[0]);

	for (int k = 0; k < n; k++) { // loop on all bubble sizes
		dBubble = dBubble_list[k];
		fprintf(stderr,"### ### ### \n");
    fprintf(stderr, "iteration %04d/%04d, dBubble[mm] %+06.5e \n",k,n,dBubble);
		fprintf(stderr,"### ### ### \n");
/* Determine Bond and Galillei number for the simulation.*/

		//double Ga = sqrt(gCons*dBubble*dBubble*dBubble*1e-9)/nuOut;
  	//double Bo = rhoOut*gCons*dBubble*dBubble*1e-6/sigma;

  	fprintf(stderr, "dBubble = %+06.3e mm\n", dBubble);
  	fprintf(stderr, "Ga      = %+06.3e\n", Ga);
  	fprintf(stderr, "Bo      = %+06.3e\n", Bo);

		rho1 = 1. [0];
  	rho2 = rho1/RHOR;
  	mu1 = 1./Ga;
  	mu2 = 1./(MUR*Ga);
  	f.sigma = 1./Bo;

  	// density
		G.x = -1.; // acceleration of gravity is set to unity... see examples/bubble.c

		run();
  }
}


event init (t = 0)
{
// Mask upper part of the domain
// that corresponds to the zone that is radially further
// from the system's axis of symmetry.
	mask (y > L0/2 ? top : none);

#if TREE
  refine (bubble_phi_axi_bary_unit (x, y, 0., Bo) > -0.1 &&
          level < MAXLEVEL);
#endif

fraction (f, bubble_phi_axi_bary_unit (x, y, 0., Bo));

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
	remove_particles(flow,  interpolate(f, x, y, z) < 0.9999);
	remove_particles(heavy, interpolate(f, x, y, z) < 0.9999);
/**
And I remove as well particles that are _too far_ from the x axis
so to reduce computational time of post-processing...*/
	remove_particles(flow,  y>1.0);
	remove_particles(heavy,  y>1.0);

/**
And I assigh how _heavy_ my particles are.*/
  foreach_particle_in(heavy)
    p().u2.z = 0.001;//0.001;  //Set relaxation timescale
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u},(double[]) {0.001,0.001,0.001,0.001},
		 maxlevel = MAXLEVEL);
}
#endif

/**
Post-processing. */

event movie(t += stepTime; t <= endTime){
  char name[256];

  view(fov=0, tx = 0, ty = -0.5, psi=-pi/2., width = 1200, height = 1200);
  draw_vof("f", lc = {0, 0, 1}, lw = 2);
  //squares("f", linear = true);
  cells(lc = {0.75, 0.75, 0.75});
  scatter(flow, s = 3);

  // mirror part of the plot, for ease of view...
  mirror(n={0,1,0}){
    draw_vof("f", lc = {0, 0, 1}, lw = 2);
    //vectors("u", scale = 0.005);
    scatter(heavy, s = 3, pc = {1, 0, 0});
  };

  snprintf(name, sizeof(name), "movie-dBubble_%+06.3e.mp4", dBubble);
  save(name);
}

event interface_tracking (t += stepTime; t <= endTime) {
  char name[256];
  snprintf(name, sizeof(name),
           "interface_location-dBubble_%+06.3e-t_%06.3f.dat",
           dBubble, t);
  FILE * fp = fopen(name, "w");
  output_facets(f, fp);
  fclose(fp);
}

event particle_tracking (t += stepTime; t <= endTime) {
  char name[256];

  // inertialess / flow particles
  snprintf(name, sizeof(name),
           "particles_flow-dBubble_%+06.3e-t_%06.3f.dat",
           dBubble, t);
  FILE * fp_flow = fopen(name, "w");

  fprintf(fp_flow,
          "# %14s %14s %14s %14s %14s %10s\n",
          "t", "x", "y", "ux", "uy", "tag");

  foreach_particle_in(flow)
    fprintf(fp_flow,
            "%+06.5e %+06.5e %+06.5e %+06.5e %+06.5e %010lu\n",
            t, x, y, p().u.x, p().u.y, p().tag);

  fclose(fp_flow);

  // inertial / heavy particles
  snprintf(name, sizeof(name),
           "particles_heavy-dBubble_%+06.3e-t_%06.3f.dat",
           dBubble, t);
  FILE * fp_heavy = fopen(name, "w");

  fprintf(fp_heavy,
          "# %14s %14s %14s %14s %14s %10s\n",
          "t", "x", "y", "ux", "uy", "tag");

  foreach_particle_in(heavy)
    fprintf(fp_heavy,
            "%+06.5e %+06.5e %+06.5e %+06.5e %+06.5e %010lu\n",
            t, x, y, p().u.x, p().u.y, p().tag);

  fclose(fp_heavy);
}