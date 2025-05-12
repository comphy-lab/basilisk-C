/**
  Single cylinder between parallel plates [Dvinsky Popel 1987](https://doi.org/10.1016/0045-7930(87)90031-4).
  Identical to [cylinder_between_plates.c], but now considering a particle with imposed angular velocity.
*/
double RADIUS = 0.15; // cylinder radius
double yShift = 0.0;
double xShift = 0.0;
double HEIGHT = 1.0+1e-8; // channel height

double nu     = 1.0; // fluid viscosity;

double uFlow     = 1.;

double TS_Omg = 25.; // TimeStepper intensity
double TS_Vel = 0.1;

coord pTor   = {0.}; // particle torque
coord pFor   = {0.}; // force
coord pTorM1 = {0.}; // torque...at previous timestep (needed for Adams Bashforth timestepper)
coord pForM1 = {0.}; // force...at previous timestep (needed for Adams Bashforth timestepper)
coord pVel   = {0.}; // translational velocity
coord pOmg   = {0.}; // angular velocity

#include "grid/quadtree.h"
#include "ghigo/src/myembed.h"
#include "fpicella/src/compute_embed_color_force_torque_RBM.h"
#include "ghigo/src/mycentered.h"
#include "fpicella/src/periodic-shift-treatment.h"

//#define FLOW (2*y*uFlow)
#define FLOW (-6.*y*y+1.5) // POISEUILLE


/**
We also define the shape of the domain. */
//#define circle (sq ((x)) + sq ((y-yShift)) - sq (RADIUS))
#define circle (sq (PS(x,xShift)) + sq (PS(y,yShift)) - sq (RADIUS))
#define channel difference(-y+HEIGHT/2.,-y-HEIGHT/2.)
#define SOLID  solid (cs, fs, intersection(channel,circle))
#define PARTICLE solid(csParticle,fsParticle,circle)


/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

int lmin = 4;  // Min mesh refinement level (l=6 is 2pt/d)
int lmax = 7; // Max mesh refinement level (l=10 is 36pt/d, l=13 is 292pt/d)
#define cmax (1.e-3) // Absolute refinement criteria for the velocity field

FILE * output_file; // global file pointer

int main ()
{
	display_control(TS_Omg,-1000,1000);
	display_control(TS_Vel,-1000,1000);
  display_control(RADIUS,0.00,0.5);
  display_control(yShift,-0.5,+0.5);
  display_control(xShift,-5.,+5.);
  display_control(lmin,1,9);
  display_control(lmax,1,9);
	display_control(uFlow,-10,10);
  /**

  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (1. [0]);

  L0 = 4.;

  X0 = Y0 = -L0/2.;
//#ifndef FLOW
//  periodic(left);
//#endif
//  periodic(top);

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

// // //	for(lmax = 8; lmax>=7; lmax -=1){
// // //		for(yShift = 0.0; yShift<= 0.2125; yShift += 0.0125)
// // //		{
// // //			RADIUS = 0.25;
// // //  		run();
// // //		}
// // //
// // //		for(RADIUS = 0.1; RADIUS <= 0.4; RADIUS += 0.05)
// // //		{
// // //			yShift = 0.;
// // //  		run();
// // //		}
// // //	}
	run();	

	fclose(output_file);
}


event init(i=0)
{
#if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */

  astats ss;
  int ic = 0;
  do {
    ic++;
		SOLID;
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
                        maxlevel = (lmax), minlevel = (lmin));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  /**
  We set the viscosity field in the event *properties*. */

  mu = muv;
	
	/**
	Set boundary conditions */
  
  u.n[left]  = dirichlet(FLOW);
  u.t[left]  = dirichlet(0.);
    p[left]  = neumann(0.);
  u.n[right] = dirichlet(FLOW);
  u.t[right] = dirichlet(0.);
    p[right] = neumann(0.);
}



event properties (i++) // refresh particle's BC on embed at each iteration...!
{
	u.n[embed]  = fabs(y)> 0.49 ? dirichlet(FLOW) : dirichlet(pVel.x -PS(y,yShift)*pOmg.x);
	u.t[embed]  = fabs(y)> 0.49 ? dirichlet(0.0 ) : dirichlet(pVel.y +PS(x,xShift)*pOmg.x);
  //required to be defined ALSO on uf (see mypoisson.h line 527...)
	uf.n[embed] = fabs(y)> 0.49 ? dirichlet(FLOW) : dirichlet(pVel.x -PS(y,yShift)*pOmg.x);
	uf.t[embed] = fabs(y)> 0.49 ? dirichlet(0.0 ) : dirichlet(pVel.y +PS(x,xShift)*pOmg.x); // pOmg.x = pOmg.y

  foreach_face()
    muv.x[] = (nu)*fs.x[];
  boundary ((scalar *) {muv});

}

event adapt (i++) // and not classic adapt, so to keep mesh identical during subiterations...
{
	SOLID;
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (lmin));
}

//event compute_particle_forces_torques (i++)
event advection_term (i++) // have this step just before starting with fluid computation
{
	// Keep forces computed from previous step (Adams Bashforth timestepper...)
	foreach()
		pForM1.x = pFor.x;
	foreach()
		pTorM1.x = pTor.x;
  //
	scalar csParticle[];
	face vector fsParticle[];
	PARTICLE;
	// allocate some auxiliary variables
	coord Fp, Fmu;
  coord Tp, Tmu;
  coord center = {0.,yShift,0.};
  //coord Omega = {oParticle,oParticle,0.}; // for the moment, it works ONLY in 2D!
  embed_color_force_RBM  (p, u, mu, csParticle, center,pOmg, RADIUS, &Fp, &Fmu);
  embed_color_torque_RBM (p, u, mu, csParticle, center,pOmg, RADIUS, &Tp, &Tmu);
  double FORCEx = Fp.x + Fmu.x;
  double FORCEy = Fp.y + Fmu.y;
  double TORQUE = Tp.x + Tmu.x;
  //fprintf(output_file,"%+6.5e %04d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
  //              L0, lmax, RADIUS, pVel.x, pOmg.x, yShift,FORCEx,FORCEy,TORQUE);
  //fflush(output_file);
// simple assignment, ok for a single particle	
	pFor.x = FORCEx;
	pFor.y = FORCEy;
	pTor.x = TORQUE;
	pTor.y = TORQUE;
	event("compute_particle_velocities");
  fprintf(stderr,"TOTO %+3.2e %02d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
                L0, lmax, RADIUS, yShift, pVel.x, pVel.y, pOmg.x,FORCEx,FORCEy,TORQUE);
}

event compute_particle_velocities (i=-1)
{
//	if(i<1){ // Explicit euler, eeded for initialization...
		foreach_dimension(){
			pVel.x += pFor.x*dt*TS_Vel;
			pOmg.x += pTor.x*dt*TS_Omg;
		}
//	}
//	else // 2-steps Adams Bashforth
//	{
//		foreach_dimension(){
//			pOmg.x += (3./2.*pTor.x-1./2.*pTorM1.x)*dt*TS_Omg;
//			pVel.x += (3./2.*pFor.x-1./2.*pForM1.x)*dt*TS_Vel;
//		}
//	}
}

/**
We look for a stationary solution. */
scalar un[]; // field that will contain previous solution...
event logfile (t += 0.01; i <= 1000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-9){ // since I'm looking for a steady solution, this must be quite small...
                           // making the simulation quite expensive
  	fprintf(output_file,"%+3.2e %02d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
                L0, lmax, RADIUS, yShift, pVel.x, pVel.y, pOmg.x,pFor.x,pFor.y,pTor.x);
  	fflush(output_file);
    return 1; /* stop */
  }
}


/**
# Force-Torque free cylinder
~~~gnuplot Free cylinder in plane Poiseuille, y=0
set grid
set xlabel "Radius"
set ylabel "U_x"
set title "Free cylinder in plane Poiseuille, y=0"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot '../Dvinsky_Popel_1987_fig_XX.dat' using ( $1):( $2) \
     with lines lw 5 lc rgb "red" title "Dvinsky Popel 1987, fig 5a",\
     'Free_Particle.dat' using ( $4==0 && $2==7 ? $3 : 1/0 ):( $4==0  && $2==7 ? +$5 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "blue" title "lmax 7",\
     'Free_Particle.dat' using ( $4==0 && $2==8 ? $3 : 1/0 ):( $4==0  && $2==8 ? +$5 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "green" title "lmax 8"
~~~

~~~gnuplot Free cylinder in plane Poiseuille, R=0.25
set grid
set xlabel "y"
set ylabel "U_x"
set title "Free cylinder in plane Poiseuille, R=0.25"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot '../Dvinsky_Popel_1987_fig_05a_R_025.dat' using ($1-0.5):( $2) \
     with lines lw 5 lc rgb "red" title "Dvinsky Popel 1987, fig 10a",\
     'Free_Particle.dat' using ( $3==0.25 && $2==7 ? $4 : 1/0 ):( $3==0.25  && $2==7 ? +$5 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "blue" title "lmax 7",\
     'Free_Particle.dat' using ( $3==0.25 && $2==8 ? $4 : 1/0 ):( $3==0.25  && $2==8 ? +$5 : 1/0 ) \
     with points pt 7 ps 1 lc rgb "green" title "lmax 8"
~~~
For the torque, it clearly is ok now, I could try to improve the computational domain
- mesh refinement
- domain size

~~~gnuplot Free cylinder in plane Poiseuille, R=0.25
set grid
set xlabel "y"
set ylabel "Omega"
set title "Free cylinder in plane Poiseuille, R=0.25"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot '../Dvinsky_Popel_1987_fig_05b_R_025.dat' using ($1-0.5):( $2) \
     with lines lw 5 lc rgb "red" title "Dvinsky Popel 1987, fig 10b",\
     'Free_Particle.dat' using ( $3==0.25 && $2==7 ? $4 : 1/0 ):( $3==0.25  && $2==7 ? +$7 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "blue" title "lmax 7",\
     'Free_Particle.dat' using ( $3==0.25 && $2==8 ? $4 : 1/0 ):( $3==0.25  && $2==8 ? +$7 : 1/0 ) \
     with points pt 7 ps 1 lc rgb "green" title "lmax 8"
~~~
For the torque, it clearly is ok now, I could try to improve the computational domain
- mesh refinement
- domain size

but honestly, it does not seem to work any better than this.
*/


