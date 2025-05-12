/*
# Moving myembed, multicolor
Using Antoon's particle structure
*/

#include "grid/quadtree.h"
#include "ghigo/src/myembed.h"
#include "fpicella/src/compute_embed_color_force_torque_RBM.h"
#include "ghigo/src/mycentered.h"
#include "fpicella/src/periodic-shift-treatment.h"

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


double xPar = 0.;
double yPar = 0.;
double omega= 0.;
double uPar = 0.;
double vPar = 0.; 
double psi  = 0.; // orientation angle

double B1 = 0.1; // 1st squirming mode intensity
double beta = 0.;// B2 = B1*beta;

#define NPARTICLES 1 // number of particles I want to work with
#include "fpicella/src/driver-myembed-particles.h"
#define CONTAINER -sq(x)-sq(y)+sq(L0/2*0.9)

int main ()
{
	display_control(xPar,-10,10);
	display_control(yPar,-10,10);
	display_control(uPar,-10,10);
	display_control(vPar,-10,10);
	display_control(omega,-10,10);
	display_control(psi,-10,10);
	display_control(B1,-10,10);
	display_control(beta,-10,10);
  /**

  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */
  
  size (1. [0]);

  L0 = 16.;

  X0 = Y0 = -L0/2.;

  stokes = true;
  DT = 1e0;//2e-5 [0];

  /**
  We initialize the grid. */

  N = 1 << (lmax);//initialize at maximum refinement, so to have _child_

  init_grid (N);
	
	mu = muv;
	
	//periodic(top);
	//periodic(right);

  run();
}
/*
### Boundary conditions
Must be declared from here, like this...
*/
event init(i = 0){
	// N.B. .n is the x direction, .t is the y direction...
	u.n[embed]   = dirichlet(surface_velocity(point,x,y,z,1));
	u.t[embed]   = dirichlet(surface_velocity(point,x,y,z,2));
	uf.n[embed]  = dirichlet(surface_velocity(point,x,y,z,1));
	uf.t[embed]  = dirichlet(surface_velocity(point,x,y,z,2));
	
	u.n[top] = dirichlet(0.);
	u.t[top] = dirichlet(0.);
 uf.n[top] = dirichlet(0.);
 uf.t[top] = dirichlet(0.);

	u.n[bottom] = dirichlet(0.);
	u.t[bottom] = dirichlet(0.);
 uf.n[bottom] = dirichlet(0.);
 uf.t[bottom] = dirichlet(0.);

	u.n[right] = dirichlet(0.);
	u.t[right] = dirichlet(0.);
 uf.n[right] = dirichlet(0.);
 uf.t[right] = dirichlet(0.);

	u.n[left] = dirichlet(0.);
	u.t[left] = dirichlet(0.);
 uf.n[left] = dirichlet(0.);
 uf.t[left] = dirichlet(0.);
}

event properties (i++)
{
	event("compute_fractions");
/*
	Adding a cylindrical container
*/
	scalar csContainer[];
	face vector fsContainer[];
	solid(csContainer,fsContainer,CONTAINER);
	foreach()
		cs[]*=csContainer[];
	foreach_face()
		fs.x[]*=fsContainer.x[];
/*
	Computing forces and adapting velocity to be force free
*/
	
	event("compute_hydro_forces");
	event("impose_force_free");

  foreach_face()
    muv.x[] = (1.0)*fs.x[];
  boundary ((scalar *) {muv});

}

event logfile (t+=0.1; i<=10000.)
{
}

//event update_particle(i++){
//	foreach_particle(){
////		p().x = xPar;
////		p().y = yPar;
////		p().u.x=uPar;
////		p().u.y=vPar;
////		p().omega.x = omega;
////		p().omega.y = omega;
////		p().omega.z = omega;
////		p().theta.x = psi;
////		p().theta.y = psi;
////		p().theta.z = psi;
//	}
//}


/*
	For visualisation purposes only
*/
event set_innner_cell_to_zero(i++){
	foreach(){
		if(cs[]==0){
			u.x[] = 0.;
			u.y[] = 0.;
		}
	}
}
