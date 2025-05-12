/**
# Passive particle in uniform flow.
identical to particle_in_uniform_flow.c

*/

//#include "embed.h"
//#include "navier-stokes/centered.h"

#include "ghigo/src/myembed.h"
#include "ghigo/src/mycentered.h"
#include "view.h"

#define RADIUS 0.25
#define SEDIMENTATION 25.
#define NPARTICLES 1
/*
Turned off, so that it is not taken into account within velocity_noslip*/
//#define imposedFlow 1
coord imposedU = {1.0,0.0,0.0}; // externally imposed flow
#include "fpicella/src/driver-myembed-particles.h"

#define MeshAdaptation 1

/*
Iterative solving...
*/
double radius = 0.;
double xShift = 0.;
double yShift = 0.;


/*
Mesh adaptation parameters */
int lmin = 3;  // Min mesh refinement level (l=6 is 2pt/d)
int lmax = 7; // Max mesh refinement level (l=10 is 36pt/d, l=13 is 292pt/d)
#define cmax (1.e-4) // Absolute refinement criteria for the velocity field

FILE * output_file; // global file pointer

int main()
{
  output_file = fopen ("Particle_Forces.dat", "w");
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (8. [0]);
  DT = HUGE [0];
  
  origin (-L0/2., -L0/2.);
//  periodic (right);
//  periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-7;
  
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

		for(radius = 0.1; radius <= 0.4; radius += 0.05)
		{
//			yShift = 0.;
  		run();
		}
//	}

	fclose(output_file);
}

scalar un[];

#define HEIGHT 1.0
#define EPS 1e-14
//#define channel difference(-x+HEIGHT/2.+EPS,-x-HEIGHT/2.-EPS) // y-aligned channel
#define channel difference(-y+HEIGHT/2.+EPS,-y-HEIGHT/2.-EPS) // x-aligned channel

event init (t = 0) {
/*Initialize particle list, structure...*/
	initialize_particles();

	foreach_particle()
		foreach_dimension()
			p().u.x = imposedU.x;
	foreach_particle()
		p().r = radius;
	foreach_particle()
		p().x = xShift;
	foreach_particle()
		p().y = yShift;
	
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
	hydro_forces_torques();
}

#if MeshAdaptation
event adapt (i++) // and not classic adapt, so to keep mesh identical during subiterations...
{
		solid(cs,fs,channel);
		compute_particle_fractions();
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (lmin));
}
#endif

/**
We check for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-4){
		foreach_particle()
			fprintf(output_file,"%+05d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e \n",lmax,p().r,p().x,p().y,p().F.x,p().F.y,p().T.x);
		fflush(output_file);
    return 1; /* stop */
	}
}

/*
Display a solution...*/
event profile (t = end) {
  printf ("\n");
  draw_vof ("cs", "fs");
  squares ("u.x", linear = false);
  vectors (u = "u", scale = 0.01, level = 8);
  save ("u.x.png");
}

/**
# Forces on the cylinder
~~~gnuplot Cylinder motion in quiescent fluid, yShift = 0
set grid
set xlabel "Radius"
set ylabel "F_x"
set title "U_x=1, y=0, quiescent fluid"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot '../../cylinder_between_plates/Dvinsky_Popel_1986_fig_02.dat' using ( $1):( $2) \
     with lines lw 5 lc rgb "red" title "Dvinsky Popel 1986, fig 2",\
     'Particle_Forces.dat' using ( $4==0 && $1==7 ? $2 : 1/0 ):( $4==0  && $1==7 ? -$5 : 1/0 ) with points pt 7 ps 2 lc rgb "blue" title "lmax 7"
~~~
*/
