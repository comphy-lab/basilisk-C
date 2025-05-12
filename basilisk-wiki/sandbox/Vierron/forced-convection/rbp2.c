/** ## Rayleigh Benard Poiseuille instability with adaptative mesh
$$
\nabla \cdot \vec{u} = 0
$$
$$
\partial_t \vec{u} + \nabla \cdot \left(\vec{u} \otimes \vec{u}
\right) = -\nabla p + Re^{-1} \nabla^{2} \vec{u} +
RaRe^{-2}Pr\theta\vec{e}_y
$$
$$
\partial_t \theta + \nabla \cdot \left( \vec{u} \theta \right) = Re^{-1}Pr^{-1} \nabla^{2} \theta
$$
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "navier-stokes/perfs.h"
#include "view.h"
#include "display.h"

/** This is the maximum and minimum level of refinement*/
int MAXLEVEL=6;
int MINLEVEL=5;

/** End time simulation */
double EndTime= 100.;

/** top limit boundary */
double ytop=1.1;

/** init vorticity */
scalar omega[];

/** init temperature field (tracer.h + diffusion.h) */
scalar T[];
scalar * tracers = {T};
mgstats mgT;

/** We need a new field to define the acceleration, viscosity and temperature diffusion. */
face vector av[], muc[], D[];

/** Physicals parameters */
double Ra, Pr, Re;

/** init domain size, grid, poisson parameter solver, periodic boundary */
int main() {
  size(4.4);
  init_grid(1<<MINLEVEL);
  DT = 0.1;
  TOLERANCE = 1e-3;
  periodic(right);
  
  Ra = 3.e3; 
  Pr = 6.4;
  Re = 2.;

  mu = muc; //viscosity
  a = av; // acceleration
  
/** When using bview we can interactively control the Reynolds, Rayleigh, Prandtl and the maximum level of refinement. */

  display_control (Re, 0.1, 1000);
  display_control(Ra, 0.1, 1e8);
  display_control(Pr,0.1,100);
  display_control(MAXLEVEL,MINLEVEL,12);
  run(); 
}

/** 
## Boundaries Conditions 
*/

T[left] = neumann(0.);

T[embed] = dirichlet(-0.5);
T[bottom] = dirichlet(0.5);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

u.n[left] = dirichlet((4.*y/ytop)*(1-(y/ytop))); // piseuille's profile

/**
## Initial conditions
Initial conditions correspond to a linear temperature
profile and a parabolic velocity profile.
*/
event init (t=0) {

/** Embedded geometry (pipe domain) */	
        vertex scalar phi[];
        foreach_vertex() 
           phi[] = ytop-y;
        boundary ({phi});
        fractions (phi, cs, fs);
        
/** Linear temperature profile + poiseuille velocity profile */        
 	foreach(){
		T[] = -((0.5+0.5)/ytop)*y+0.5;
		u.x[] = (4.*y/ytop)*(1-(y/ytop));
 	}
        boundary ({T,u});
}

/** ## Set viscosity */
event properties (i++) {
	foreach_face()
		muc.x[] = fm.x[]/Re;	
	boundary ((scalar*){muc});
}

/** ## Set acceleration */
event acceleration (i++) {
	foreach_face(y)
		av.y[] += (Ra*Pr/sq(Re))*(T[] + T[0,-1])/2.;
}

/** ## Set temperature diffusion */
event tracer_diffusion (i++) {
	foreach_face()
		D.x[] = fm.x[]/(Re*Pr);
	boundary ((scalar*){D});
	mgT = diffusion (T, dt, D);
  	boundary({T});
}

/** 
## Adaptative mesh on embedded geometry and vorticity
*/
#if 1
event adapt (i++) {
  	vorticity (u, omega);
	adapt_wavelet ({cs,omega}, (double[]){1.e-2,5.e-2}, MAXLEVEL, MINLEVEL);
	
/** refine temperature boundary layers top and bottom */
	refine(y<0.1 && level < MAXLEVEL+1);
	refine(y>ytop-0.1 && level < MAXLEVEL+1);
}
#endif

/**## Output data */
/** output_field temperature and velocity */
event output (t+=0.5){
  	char names[80];
	sprintf(names, "data%g", t);
	FILE * fout = fopen (names, "w");
	output_field({T,u}, fout);
}
	

/** 
## Output movies
*/
event movies (i++) {

/** grid evolution */
	scalar l[];
	foreach(){
		l[] = level;
	}
	output_ppm (l, file="grid.mp4", n=512, box = {{0.,0.},{L0,ytop}});
	
/** temperature field evolution */
	output_ppm (T, file="temp.mp4", n=1024, linear=true, box = {{0.,0.},{L0,ytop}});
	
#if 1
  	scalar omega[];
  	vorticity (u, omega);
	box();
	view (quat = {0.000, 0.000, 0.000, 1.000},
      	fov = 30, near = 0.01, far = 1000,
      	tx = -0.469, ty = -0.088, tz = -1.359,
      	width = 1920, height = 1080, samples = 4);
	squares("omega",linear=true);
	save("omega.mp4");

	clear();

	box();
	view (quat = {0.000, 0.000, 0.000, 1.000},
      	fov = 30, near = 0.01, far = 1000,
      	tx = -0.469, ty = -0.088, tz = -1.359,
      	width = 1920, height = 1080, samples = 4);
	squares("T", min=-0.5, max=0.5, linear=true);
	cells();
	save("temperature.mp4");
#endif
}

/** dump */
event snapshot (t=EndTime) {
	dump();
}

/**
## Results

![Temperature field.](RBP/temperature.mp4)
![Animation of the level of refinement.](RBP/grid.mp4)
*/


