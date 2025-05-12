/**
# The idealized case of ventilator operation in the orchard
*/
/**
this code mostly is adapted from [Vincent](http://basilisk.fr/sandbox/vheusinkveld/idealized.h) work. The [fan](http://basilisk.fr/sandbox/YiDai/fan.h) function is changed. Other functions are intact. 
*/


/** Include required libraries */ 
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "grid/octree.h" 		// For 3D
#include "view.h"			// For bview
#include "navier-stokes/centered.h"     // Navier stokes 
#include "tracer.h"			// Tracers
#include "diffusion.h"			// Diffusion 

/** Global variables */
int minlevel, maxlevel;         	// Grid depths
double meps, eps;			// Maximum error and error in u fields
double TEND = 1140;

#include "physics.h"			// Physics of the simulation 
#include "fan.h"			// Include a fan
#include "output_slices.h"
#include "lambda2.h"
// #include "diagnostics.h"		// Perform diagnostics

/** Initialisation */
int main() {	
    minlevel = 4;
    maxlevel = 9;

    L0 = 700.;
    X0 = Y0 = Z0 = 0.;

    // Possibility to run for variable changes
    // for(double tempVar=90; tempVar<110; tempVar+=3) {
    // //for(rot.theta=90*M_PI/180; rot.theta<106*M_PI/180.; rot.theta+=3*M_PI/180) {
	
    //     rot.theta = tempVar*M_PI/180;

    init_grid(1<<6);
    a = av; 

    foreach_dimension() {
        u.x.refine = refine_linear;  		// Momentum conserved 
    }

	fan.prolongation = fraction_refine;		// Fan is a volume fraction
	p.refine = p.prolongation = refine_linear;
	b.gradient = minmod2; 				// Flux limiter 

    rot.phit = 2*M_PI/60;
    // rot.ST_phit = M_PI;
    rot.ST_phit = -90.; // full rotation

  	meps = 10.;					// Maximum adaptivity criterion
	DT = 10E-5;					// For poisson solver 
    TOLERANCE=10E-6;				// For poisson solver 
	CFL = 0.8;					// CFL condition

    run();						// Start simulation 

}


/** Initialisation */
event init(t=0) {
    rot.fan = true;		// Yes we want a fan
    rot.rotate = true;		// If we want it to rotate 
    rot.start = 60;
    rot.stop = 1020;
    eps = .5;
    
    init_physics();

    if(rot.fan) {
        init_rotor();
	    rotor_coord();
    }
        
    while(adapt_wavelet((scalar *){u,b},(double []){eps,eps,eps,0.35*9.81/273},maxlevel,minlevel).nf) {
	foreach() {
	    b[] = STRAT(y);
        u.x[] = WIND(y);
	}
	rotor_coord();
    }
}

/** Return to standard tolerances and DTs for poisson solver */ 
event init_change(i=10) {
    TOLERANCE=10E-3;
    DT = .5;
}

/** Adaptivity */
event adapt(i++) {
    adapt_wavelet((scalar *){fan,u,b},(double []){0.01,eps,eps,eps,.35*9.81/273},maxlevel,minlevel);
}

/** Progress event */
event progress(t+=5) {
    fprintf(stderr, "i=%d t=%g p=%d u=%d b=%d \n", i, t, mgp.i, mgu.i, mgb.i);
}
event mov (t += 0.5) {
#if (dimension == 2)
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  // draw_vof ("cs", "fs", filled = -1, fc = {0,0,0});
  // draw_vof ("cs", "fs");
  squares ("omega", linear = true, map = cool_warm);
  mirror ({0,-1})
  cells();
#elif (dimension == 3)
  scalar l2[];
  lambda2 (u, l2);
  view (fov = 30, theta = 0.5, phi = 0.4, 
	tx = -0.1, ty = 0.1, bg = {65./256,157./256,217./256},
	width = 1080, height = 1080);
  box();
  isosurface ("l2", -0.01);
  translate (y = -5){
    squares ("u.x", n = {0,1,0}, alpha = 5.,
	     min = -1.2, max = 2, map = cool_warm);}
  translate (z = -L0/2.){
    squares ("b", n = {0,0,1}, alpha = L0/2.,
	     min = -0.1, max = 0.6, map = cool_warm);}
  translate (x = -L0/2.){
    squares ("b", n = {1,0,0}, alpha = L0/2.,
	     min = -0.1, max = 0.6, map = cool_warm);}
#endif
  save ("lam_3D.mp4");
}


char dir_slices[61];

event slices(t += 1) {
    char nameSlice[91];
    coord slice = {1., 0., 1.};
    int res = L0/2;

    for(double yTemp = 0.5; yTemp<=5; yTemp+=0.5) {
	    slice.y = yTemp/L0;

    	snprintf(nameSlice, 90, "%st=%05gy=%03g", dir_slices, t, yTemp*10);
    	FILE * fpsli = fopen(nameSlice, "w");
    	output_slice(list = (scalar *){b}, fp = fpsli, n = res, linear = true, plane=slice);
    	fclose(fpsli);
    }
}


/** End the simulation */
event end(t=TEND) {
}
