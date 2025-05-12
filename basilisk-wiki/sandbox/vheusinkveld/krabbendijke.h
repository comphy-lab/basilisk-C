/** 
## Intro
This is the main.c file of the Krabbendijke case. Here all main simulation settings are put.

* The physics functions and events can be found [here](http://basilisk.fr/sandbox/vheusinkveld/afm/krabbendijke/physics.h)
* The fan functions and events can be found [here](http://basilisk.fr/sandbox/vheusinkveld/afm/krabbendijke/fan.h)
* The diagnostics functions and events can be found [here](http://basilisk.fr/sandbox/vheusinkveld/afm/krabbendijke/diagnostics.h)

For the imported centered Navier-Stokes formulation, see [this page](http://basilisk.fr/src/navier-stokes/centered.h)

## Overview of the simulation
### main.c (here named krabbendijke.h)

* Simulation run time (TEND)
* Size (L0) and resolution (2^maxlevel=max gird points in 1D) 
* For loop around run() function such that variables can be varied. 
* The simulation is initialized (also using functions from fan.h and physics.h). In this initialization the grid is adapted to the initialzed velocity and buoyancy fields. 
* At every time step adaptivity is taken care of and some information is printed.

### [physics.h](http://basilisk.fr/sandbox/vheusinkveld/afm/krabbendijke/physics.h)

* Physics init functions sets boundary conditions and velocity and buoyancy field (after which the grid is adapted to them in main.c)
* Subgrid scale closure is turned on via SGS.h
* Inversion (STRAT) and wind (WIND) speed are defined.  
* The surface temperature is determined based on the lowest 2 grid levels, based on this temperature a surface flux is calculated. This is treaded as a source in the diffusion routine for the buoyancy.
* The acceleration event applied 'gravity' via the buoyancy field
* In the outer part of the domain wind and buoyancy are forced to their initialized profiles with a relaxation time

### [fan.h](http://basilisk.fr/sandbox/vheusinkveld/afm/krabbendijke/fan.h)

* The rot structure is defined, containing properties of the fan/rotor (i.e. height, power, rotation time)
* Fan is rotated at every simulation step and corresponding properties updated (i.e. normal vector, volume fraction field)
* Forcing is applied at every simulation step.

### [diagnostics.h](http://basilisk.fr/sandbox/vheusinkveld/afm/krabbendijke/diagnostics.h)

* dia structure is defined containing diagnosed properties (i.e. kinetic energy in the system, work done by the fan)
* out structure is defined containing output information (i.e. folder names, file name prefixes)
* bvsets structure is defined containing settings for bview, the movie output
* A 'case file' is written away containing all the main properties of a given run
* Simple diagnostics are outputted frequently (user defined) 
* Line windmeasurments are done every second at 3 m height
* Line temperature (buoyancy) measurments are done every second at 1 m height for a given angle wrt wind direction
* Temperature slices are outputted every 5 seconds
* Movies are generated in which eddies are visualized and gird depth, buoyancy and wind velocity are projected on the background.

## The code of main.c or krabbendijke.h

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
double meps, eps;			// Maximum and normal refinement criteria
double TEND = 1560.;

char sim_ID[] = "krab";		        // Simulation identifier
char sim_var[] = "angle";  		// Notes if a variable is varied over runs

#include "physics.h"			// Physics of the simulation 
#include "fan.h"			// Include a fan
#include "diagnostics.h"		// Perform diagnostics

/** Initialisation */
int main() {	
    minlevel = 4;
    maxlevel = 9;

    L0 = 700.;
    X0 = Y0 = Z0 = 0.;

    // Possibility to run for variable changes
    for(double tempVar=300; tempVar<301; tempVar+=100) {
	
	//rot.theta = tempVar*M_PI/180;
	rot.phit = 2*M_PI/tempVar;                      // Hence a 300 second rotation time

        init_grid(1<<6);
	a = av; 

        foreach_dimension() {
	    u.x.refine = refine_linear;  		// Momentum conserved 
	}

	fan.prolongation = fraction_refine;		// Fan is a volume fraction
	p.refine = p.prolongation = refine_linear;
	b.gradient = minmod2; 				// Flux limiter 

  	meps = 10.;					// Maximum adaptivity criterion
	DT = 10E-5;					// For poisson solver 
        TOLERANCE=10E-6;				// For poisson solver 
	CFL = 0.8;					// CFL condition

	sim_dir_create();				// Create relevant dir's
	out.sim_i++;					// Simulation iteration
 
    	run();						// Start simulation 

    }
}

/** Initialisation */
event init(t=0) {
    rot.fan = true;		// Yes we want a fan
    rot.rotate = true;		// If we want it to rotate 
    rot.start = 30.;            // start time of rotor forcing
    rot.stop = 1260.;           // stop time of rotor forcing

    rot.phi = 0;	        // Reset for repeated different runs
    eps = .5;
    
    init_physics();             // See imported phyisics module for this function

    if(rot.fan) {               // See imported fan module for these functions
        init_rotor();
	rotor_coord();
    }
    // Adapt grid to initial buoyancy and velocity fields.    
    while(adapt_wavelet((scalar *){u,b},(double []){eps,eps,eps,0.25*9.81/273},maxlevel,minlevel).nf) {
	foreach() {
	    b[] = STRAT(y);  // physics.h
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
    adapt_wavelet((scalar *){fan, u,b},(double []){0.01,eps,eps,eps,.25*9.81/273},maxlevel,minlevel);
}

/** Progress event */
event progress(t+=5) {
    fprintf(stderr, "i=%d t=%g p=%d u=%d b=%d \n", i, t, mgp.i, mgu.i, mgb.i);
}

event dumpfields(t=120; t+=120) {
    char nameDump[90];
    snprintf(nameDump, 90, "./%s/fielddump", out.dir);
    dump(file = nameDump, list = all);
}

/** End the simulation */
event end(t=TEND) {
}
