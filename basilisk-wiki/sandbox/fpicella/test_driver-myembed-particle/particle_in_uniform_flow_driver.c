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
#define imposedFlow 1
coord imposedU = {1.0,1.0,0.0}; // externally imposed flow
#include "fpicella/src/driver-myembed-particles.h"




int main()
{
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (2. [0]);
  DT = HUGE [0];
  
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-7;
  
  //for (N = 16; N <= 64; N *= 2)
	N = 32;
    run();
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
	
  /**
  Viscosity is unity. */
  
  mu = fm;

  /**
  The channel geometry is defined using Constructive Solid Geometry. */  

	solid(cs,fs,channel);
 /**
	Add as well the presence of the particle. */
	compute_particle_fractions();

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

/**
We check for a stationary solution. */

event logfile (t += 0.1; i <= 1000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
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
