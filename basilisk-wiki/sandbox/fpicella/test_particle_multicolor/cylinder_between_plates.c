/**
	# Single cylinder between parallel plates
	Dvinsky Popel 1987 is the goal.
	adapted from Ghigo's myembed.h
*/
#define p_n (3)
/**
Quantities imposed on the particle
*/
double uP=1.;
double vP=0.;
double wP=0.; // imposed angular velocity of the particle...
double yShift = 0.15;
double RADIUS = 0.25;

#include "grid/quadtree.h"
#include "ghigo/src/myembed.h"
#include "fpicella/src/compute_embed_color_force_torque_RBM.h"
#include "ghigo/src/mycentered.h"
#include "ghigo/src/myembed-moving-multicolor.h"

//#include "view.h"

// External flow
//#define FLOW (0.)
//#define FLOW (-6.*y*y+1.5) // POISEUILLE
//#define FLOW (2.*y+0.) // COUETTE

/**
## Reference solution */
#define d    (2.0*RADIUS)    // Diameter of the cylinder
#define h    (1.00001)  // Width of the channel

#define nu   (1.0)

#define uref 1.
#define tref 1.

/**
We also define the shape of the domain. The first two particles are in
fact the channel, while the other particle is the cylinder. */

#define circle(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))
#define wall(y,w)   ((y) - (w)) // + over, - under

// Channel
double p0_phi (double xc, double yc) // height of the channel...
{
  //return intersection (-(wall (yc,  (h)/2.)), (wall (yc, -(h)/2.)));
  return intersection (-(wall (yc,  (h)/2.)), 1.);
}
double p1_phi (double xc, double yc) // length of the channel...
{
  //return intersection (-(wall (xc, (l))), (wall (xc, (-l))));
  return intersection (1., (wall (yc, -(h)/2.)));
}
// Cylinder
double p2_phi (double xc, double yc)
{
  return circle (xc, yc);
}

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

int lmin = 4;  // Min mesh refinement level (l=6 is 2pt/d)
int lmax = 8; // Max mesh refinement level (l=10 is 36pt/d, l=13 is 292pt/d)
#define cmax (1.e-6*(uref)) // Absolute refinement criteria for the velocity field

FILE * output_file; // global file pointer

int main ()
{
	display_control(RADIUS,0.00,0.5);
	display_control(yShift,0.00,0.5);
	display_control(lmin,1,9);
	display_control(lmax,1,9);
  /**

  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */
  
  size (1. [0]);
  
  L0 = 8.;
  
  X0 = Y0 = -L0/2.;
  
  /**
  We turn off the advection term. The choice of the maximum timestep
  and of the tolerance on the Poisson and viscous solves is not
  trivial. This was adjusted by trial and error to minimize (possibly)
  splitting errors and optimize convergence speed. */
  
  stokes = true;
  DT = 1e-2;//2e-5 [0];

  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6*(uref);
  
  /**
  We initialize the grid. */

#if !TREE
  N = 1 << (lmax);
#else
  N = 1 << (lmin);
#endif // !TREE
  init_grid (N);
    
  output_file = fopen ("Particle.dat", "w");

//	for(lmax = 7; lmax>=7; lmax -=1){
//		for(yShift = 0.0; yShift<= 0.2125; yShift += 0.0125)
//		{
//			RADIUS = 0.25;
//  		run();
//		}
//
//		for(RADIUS = 0.1; RADIUS <= 0.4; RADIUS += 0.05)
//		{
//			yShift = 0.;
//  		run();
//		}
//	}
//
//	fclose(output_file);

	run();
}

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = (nu)*fm.x[];
  boundary ((scalar *) {muv});
}

/**
## Initial conditions */

event init (i = 0)
{
  /**
  We set the viscosity field in the event *properties*. */

  mu = muv;

  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */
    
#if ORDER2
  for (scalar s in {u, p, pf})
    s.third = false;
#else
  for (scalar s in {u, p, pf})
    s.third = true;
#endif // ORDER2

  /**
  We use a slope-limiter to reduce the errors made in small-cells. */
  
#if SLOPELIMITER
  for (scalar s in {u, p, pf}) {
    s.gradient = minmod2;
  }
#endif // SLOPELIMITER
    
  /**
  We initialize the embedded boundary.

  We define the shape of the three particles. */

  pl[0].phi = p0_phi;
  pl[1].phi = p1_phi;
  pl[2].phi = p2_phi;

  /**
  Next we define the particles' physical parameters. The channel is
  defined using the default parameters. */

//  // Ratio of solid and fluid density
//  pl[2].r = (rr);
//  // Particle volume
//  pl[2].v = (1.);
//  // Particle moment of inertia
//  foreach_dimension()
//    //pl[2].i.x = (p_moment_inertia_cylinder ((d), pl[2].r));
//    pl[2].i.x = (rr);

  /**
  We then define the particles' initial position. */

  pl[2].c.x = 0.;
  pl[2].c.y = yShift;

  /**
  We then define the particles' imposed velocity. */

  pl[2].u.x = uP;
  pl[2].u.y = vP;
  pl[2].w.x = wP;
  pl[2].w.y = wP; // must impose on both components in 2D, sorry!
	
	///**
	//Define initial velocity */
	//foreach()
	//	u.x[] = FLOW;
  
#if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */
  
  astats ss;
  int ic = 0;
  do {
    ic++;
    p_shape (cs, fs, pl);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
                        maxlevel = (lmax), minlevel = (lmin));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  
  p_shape (cs, fs, pl);
#ifdef FLOW
  u.n[left] = dirichlet(FLOW);
  u.t[left] = dirichlet(0.);
  p[left] = neumann(0.);
  u.n[right] = dirichlet(FLOW);
  u.t[right] = dirichlet(0.);
    p[right] = neumann(0.);
#endif

/**
We remove the Neumann pressure boundary condition which is responsible for instabilities.
*/

//  for (scalar s in {p, pf, g}) {
//    s.neumann_zero = true; // turned off
//   // s.neumann_zero = false; // turned on
//  }

}

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++) 
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (lmin));

  /**
  We do not need here to reset the embedded fractions to avoid
  interpolation errors on the geometry as the is already done when
  moving the embedded boundaries. It might be necessary to do this
  however if surface forces are computed around the embedded
  boundaries. */
}
#endif // TREE

/**
## Outputs */
event particle_location(i++){
  /**
  We then define the particles' initial position. */

  pl[2].c.x = 0.;
  pl[2].c.y = yShift;
}

event particle_force_torque (i++) 
{
  coord Fp, Fmu;
  coord Tp, Tmu;
  embed_color_force_RBM  (p, u, mu, pl[2].col, pl[2].c, pl[2].w, RADIUS, &Fp, &Fmu);
  embed_color_torque_RBM (p, u, mu, pl[2].col, pl[2].c, pl[2].w, RADIUS, &Tp, &Tmu);
	/**
	Store particles' force and torque within the acceleration variables
	*/
  pl[2].au.x = Fp.x + Fmu.x;
  pl[2].au.y = Fp.y + Fmu.y;
  pl[2].aw.x = Tp.x + Tmu.x;
  fprintf(stderr,"TOTO %+3.2e %02d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
                L0, lmax, RADIUS,pl[2].c.y,
								pl[2].u.x,  pl[2].u.y,  pl[2].w.x,
								pl[2].au.x, pl[2].au.y, pl[2].aw.x
				 );
}

scalar un[]; // field that will contain previous solution...
event logfile (i++; i<500) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-5){
  	fprintf(output_file,"%+3.2e %02d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
                L0, lmax, RADIUS,pl[2].c.y,
								pl[2].u.x,  pl[2].u.y,  pl[2].w.x,
								pl[2].au.x, pl[2].au.y, pl[2].aw.x
		);
		fflush(output_file);
    return 1; /* stop */
  }
}


/**
# Forces on the cylinder

~~~gnuplot Cylinder motion in quiescent fluid, RADIUS = 0.25, variable yShift, TORQUE
set grid
set xlabel "y"
set ylabel "T"
set title "U_x=1, y=0, quiescent fluid, R=0.25"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot '../Dvinsky_Popel_1986_fig_04_R_025.dat' using ($1-0.5):( $2) \
     with lines lw 5 lc rgb "red" title "Dvinsky Popel 1986, fig 4",\
     'Particle.dat' using ( $3==0.25 && $2==7 ? $6 : 1/0 ):( $3==0.25 && $2==7 ? -$9 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "blue" title "lmax 7",\
     'Particle.dat' using ( $3==0.25 && $2==8 ? $6 : 1/0 ):( $3==0.25 && $2==8 ? -$9 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "black" title "lmax 8",\
     'Particle.dat' using ( $3==0.25 && $2==9 ? $6 : 1/0 ):( $3==0.25 && $2==9 ? -$9 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "green" title "lmax 9"\
~~~
For the torque, it clearly suffers for strong confinement.
I should _definitely_ try for higher resolution.
*/
