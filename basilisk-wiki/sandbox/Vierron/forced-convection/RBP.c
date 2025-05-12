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

/** ## Rayleigh Benard instability with adaptative mesh */
#include "boussinesq_poiseuille.h"
#include "navier-stokes/perfs.h"
#include "view.h"
#include "display.h"


//#define MINLEVEL 6
//#define MAXLEVEL 7

int MAXLEVEL=7;
int MINLEVEL=5;
double EndTime= 100.;
double ytop=2.;

scalar omega[];

int main() {
  size(8);
  init_grid(1<<MINLEVEL);
  Ra = 3.e3; 
  Pr = 0.72;
  Re = 0.5;
  DT = 0.1;
  TOLERANCE = 1e-3;
  periodic(right);
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
//T[right] = neumann(0.);

T[embed] = dirichlet(-0.5);
T[bottom] = dirichlet(0.5);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

u.n[left] = dirichlet((4.*y/ytop)*(1-(y/ytop)));
p[left] = neumann(0.);
pf[left] = neumann(0.);

/**
## Initial conditions
Initial conditions correspond to a linear temperature
profile and a parabolic velocity profile.
*/
event init (t=0) {
	
        vertex scalar phi[];
        foreach_vertex() 
           phi[] = ytop-y;
        boundary ({phi});
        fractions (phi, cs, fs);
        
 	foreach(){
		T[] = -((0.5+0.5)/ytop)*y+0.5;
		u.x[] = (4.*y/ytop)*(1-(y/ytop));
 	}
        boundary ({T,u});
}

/** 
## Adaptative mesh 
*/
#if 1
event adapt (i++) {
  	vorticity (u, omega);
	adapt_wavelet ({cs,omega}, (double[]){1.e-2,5.e-2}, MAXLEVEL, MINLEVEL);
	refine(y<0.1 && level < MAXLEVEL+1);
	refine(y>ytop-0.1 && level < MAXLEVEL+1);
}
#endif

event out (t+=0.5){
  	char names[80];
	sprintf(names, "data%g", t);
	FILE * fout = fopen (names, "w");
	output_field({T,u}, fout);
}
	

/** 
## Output Data t+=1.; t <= EndTime)
*/
event output (i++) {

	scalar l[];
	foreach(){
		l[] = level;
	}
	output_ppm (l, file="grid.mp4", n=512, box = {{0.,0.},{L0,ytop}});
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

/** 
## Calcul profil de temperature + velocity
*/
event pro (t=EndTime) {
	//profile({T,u},"profils",0, ytop);
	dump();
	return 0;
}

/**
## Results

![Temperature field.](RBP/temperature.mp4)
![Animation of the level of refinement.](RBP/grid.mp4)
*/


