/** 
# Intro 

This is a minimal example for implementing a 'fan' after geometry examples that can be found [here](http://basilisk.fr/sandbox/Antoonvh/puck.c). The mentioned function is updated to be applicable for general phi and theta. Also this script adds a velocity based on a stated fan power. All functions can be found at the end of this file.

# Power to velocity forcing and expected exit velocity
The implementation is based on a momentum forcing. An amount of work per time is done on the flow with a certain power, P. Using the work-energy theorem, which states that the work done by all acting forces on a particle is equal the change in kinetic energy, this can be expressed as $\Delta E_{kin}=P\Delta t$. This directly implies that,

$$\vec{v}^2_{t+\Delta t} - \vec{v}^2_{t} = \frac{2\textrm{P}\Delta t}{m}\hat{n}^2$$

which in turn implies for every spatial dimension,

$$ v^2_{i,t+\Delta t} - v^2_{i,t} = \frac{2\textrm{P}\Delta t}{m} n_i^2 $$

From these equations the direction of $\hat{n}$ is not yet defined and $|\hat{n}|=1$. When it is assumed that the force applied to the flow is normal to the plane in which the fan blades rotate then $\hat{n}$ is the normal vector of this plane. Rearranging and taking the square root results in,
$$	\vec{v}_{i,t+\Delta t} = \sqrt{\vec{v}^2_{i,t} + \frac{2\textrm{P}\Delta t}{m}\hat{n}^2_i } $$

To estimate the velocity to which air is accelerated by this forcing we assume that it is applied over a finite time from $t_{start}$ to $t_{end}$, moving a distance $w$ in the process:
$$		w=\int_{t_{start}}^{t_{end}} v\left(t\right)\textrm{dt}$$
Assuming a zero entrance velocity, $v$ is given at every instance as:
$$\textrm{E}_{\textrm{kin}}=\textrm{P}\Delta t=0.5mv^2$$
$$v\left(t\right)=\sqrt{2\textrm{P}t/m}$$
P, $w$ and m are set by the user thus this leaves two equations and two unknowns. Solving for the exit velocity gives:
$$	v_{end} = \left(\frac{3\textrm{P}w}{m}\right)^{1/3}$$


# The code
First we include relevant files for 2D navier stokes solving and geometric fun.
*/
//#include "grid/octree.h" 		// For 3D --> to demanding for server
#include "navier-stokes/centered.h"     // Navier stokes 
#include "fractions.h"	      		// Needed to compute fan fractions

/** Required functions and structures */
void init_rotor();
void rotor_update();
void rotor_coord();
void rotor_forcing();

struct sRotor {	
	double rampT;		// Time to start up rotor
	double P, Prho;		// Power, powerdensity 
	double R, W, A, V;	// Diameter, Thickness, Area ,Volume
	double diaVol;		// Diagnosed rotor volume 
	double x0, y0, z0;	// Origin of rotor
    	double xt, yt, zt;      // Translation of origin per rotate event
        double theta, phi;	// Polar and Azimuthal angle 
   	double thetat, phit;    // Addition to angles per rotate event
   	double Work;            // Work done by rotor
   	double cu;              // Characteristic velocity
   	bool rotate;            // Rotation yes or no
	coord nf, nr;		// Normal vector fan, rotation 
};

/** Some global vars */
int minlevel, maxlevel;         // Grid depths
double eps;			// Adaptivity criterion
struct sRotor rot;  		// Rotor details structure 
scalar fan[];			// Fan volume fraction
scalar lev[]; 			// for movie diagnostic level 
scalar Ekin[];  		// for movie kinetic energy 

/** Starting with the most loved c function */

int main() {	
	minlevel = 3;
  	maxlevel = 8;

   	L0 = 10.;
   	X0 = Y0 = Z0 = 0.;

    	init_grid(1<<7);

	foreach_dimension(){
		periodic (left);
		u.x.refine = refine_linear; 		// Momentum conserved
        }

	fan.prolongation = fraction_refine;		// Fan is a volume fraction
	p.refine = p.prolongation = refine_linear;

	DT = 10E-5;					// For poisson solver 
        TOLERANCE=10E-6;				// For poisson solver 
	CFL = 0.5;					// CFL condition

    	run();						// Start simulation 	
}

/** Initialisation */
event init(t = 0){
	rot.rotate = true; 				// Turn on rotation
	init_rotor();					// See init function for details 
	fan.prolongation = fraction_refine;		// Tell basilisk that fan is volumetric field
	refine (fan[] > 0. && level < maxlevel); 	// Refine where fan is
	eps = 0.07*rot.cu;				// Wavelet estimator based on flow velocity
}

/** Reset to defaults for poisson solver */
event reset_pois(i = 10){
  DT = 1;
  TOLERANCE = 10E-3;
}

/** Adaptivity */
event adapt(i++) {
	adapt_wavelet((scalar *){fan,u},(double []){0.,eps,eps},maxlevel,minlevel);
}

/** Forcing by the rotor */
event forcing(i = 1; i++) {
	rotor_coord();
	rotor_forcing();	
}

/** Rotate the rotor */
event rotate(t+=0.1) {
    if(rot.rotate) { 
        // Change center  
        rot.x0 += rot.xt;
        rot.y0 += rot.yt;
        rot.z0 += rot.zt;

        // Change angles 
        rot.theta += rot.thetat;
        rot.phi += rot.phit;

        rotor_update();
    }
}

/** Visual output */
event movies(t+=1){
    fprintf(stderr, "t=%g\n", t);
    foreach(){
   	lev[] = level;
 	Ekin[] = sqrt(sq(u.y[]) + sq(u.x[]));
    }
    boundary({lev, Ekin});

    output_ppm (lev, file = "ppm2mp4 level.mp4", n = 512, linear = false, min = minlevel, max = maxlevel);
    output_ppm (Ekin, file = "ppm2mp4 ekin.mp4", n = 512, linear = true);


}

/** End the simulation */
event end(t=180){
}
/** 
# Results
![Kinetic energy density](fan/ekin.mp4)
![Grid depth](fan/level.mp4)
*/

/** 
# Functions 
*/

/** Function returning the sRotor structure, includign default properties */
void init_rotor() {
	rot.Work = 0.;
    	if(!rot.rampT)
	    	rot.rampT = 1.; // fan at full thrust after 1 second, linear ramp
    	if(!rot.R)
	    	rot.R = L0/80.;     
    	if(!rot.W)
	    	rot.W = rot.R/5;    
    	if(!rot.Prho)                  
     		rot.Prho = L0/10.;		
   	if(!rot.x0)
    		rot.x0 = L0/2.;
    	if(!rot.y0)
	    	rot.y0 = L0/2.;
    	if(!rot.z0){
        #if dimension == 2
            	rot.z0 = 0.;
        #elif dimension == 3
            	rot.z0 = L0/2.;
        #endif
        }
    	if(!rot.theta)
	    	rot.theta = 90.*M_PI/180.;      // Polar angle
    	if(!rot.phi)
	    	rot.phi = 90.*M_PI/180.;        // Azimuthal angle 
  
    	if(rot.rotate) {
                // determines for every rotate event how much needs to be added to x, y, z, tehta, phi.
        	rot.xt = 0;
        	rot.yt = 0;
        	rot.zt = 0;
        	rot.thetat = 0.;
       	 	rot.phit = -0.1*M_PI/180.; // note that this rate is linked to the frequency of the rotate event
    	} else {
       		rot.xt = 0;
        	rot.yt = 0;
        	rot.zt = 0;
        	rot.thetat = 0.;
        	rot.phit = 0.;
    }
	rotor_update();
}

/** Updating relevant rotor variables */
void rotor_update() {
   	// Normal vectors to the fan plane, note that z and y are turned around for consistency regarding polar coordinates (basilisk uses y for the vertical)
    	rot.nf.x = sin(rot.theta)*cos(rot.phi);
	rot.nf.y = sin(rot.theta)*sin(rot.phi);
	rot.nf.z = cos(rot.theta);

	rot.nr.x = sin(rot.theta)*cos(rot.phi);
   	rot.nr.y = sin(rot.theta)*sin(rot.phi);
    	rot.nr.z = cos(rot.theta);

   	#if dimension == 2	
		rot.A = 2*rot.R*rot.W;
	#elif dimension == 3    	
		rot.A = sq(rot.R)*M_PI;      
	#endif

        // Volume is slice of a sphere
  	#if dimension == 2
		rot.V = 1.*rot.A;
	#elif dimension == 3 
		rot.V = 4.*M_PI*pow(rot.R,3.)/3. - 
			2*M_PI*pow(rot.R-rot.W/2., 2.)/3.*(2*rot.R + rot.W/2.);
	#endif

	rot.P = rot.V*rot.Prho;
    	rot.cu = pow(3*rot.Prho*rot.W, 1./3.);
}

/** Function returning the volume fractions of a fan object */
void rotor_coord() {
    	scalar sph[], plnu[], plnd[];
   	fraction(sph, -sq((x - rot.x0)) - sq((y - rot.y0)) - sq((z - rot.z0)) + sq(rot.R)); // sphere
  	fraction(plnu,  rot.nr.x*(x - rot.x0) + rot.nr.y*(y - rot.y0) + rot.nr.z*(z - rot.z0) + rot.W/2.); // upper plane
   	fraction(plnd, -rot.nr.x*(x - rot.x0) - rot.nr.y*(y - rot.y0) - rot.nr.z*(z - rot.z0) + rot.W/2.); // lower plane

	foreach () {
    		fan[] = sph[]*plnu[]*plnd[]; // estimate for the overlappnig volumes, not that this is not exact at the boundaries of the 'fan', this is corrected for later in the focring by comparing diagnosed to theoretical volume. 
   	}
	boundary({fan}); // update boundaries
}

/** Function returning new velocities based on a rotor forcing.
This is a function of powerdensity, width, direction, ramp-time, and diagnosed volume */
void rotor_forcing(){
	double tempW = 0.;
	double w, wsgn, damp, usgn, utemp, corrP;
	foreach(reduction(+:tempW)) {		
		if(fan[] > 0.) {
			foreach_dimension() {
			wsgn = sign(rot.nf.x*u.x[]) + (sign(rot.nf.x*u.x[]) == 0)*sign(rot.nf.x); // Get the direction of the flow relative to forcing direction 
			damp = rot.rampT > t ? t/rot.rampT : 1.;                                  // Linear ramp for starting up the forcing
			corrP = rot.diaVol > 0. ? rot.V/rot.diaVol : 1.;                         // Correction if diagnosed volume differs from theoretical volume
			w = wsgn*fan[]*damp*sq(rot.nf.x)*(2./rho[])*(corrP*rot.P/rot.V)*dt;      // additional kinetic energy
			tempW += 0.5*rho[]*w*dv(); // save work done to asses total work in diagnostics

			// New kinetic energy
			utemp = sq(u.x[]) + w;
                        // New sign of the velocity
			usgn = 	  1.*(u.x[] >= 0)*(utemp > 0) +
			    	 -1.*(u.x[] >= 0)*(utemp < 0) +
		 		  1.*(u.x[] <  0)*(utemp < 0) +
				 -1.*(u.x[] <  0)*(utemp > 0); 
                        // New velocity
			u.x[] = usgn*sqrt(fabs(utemp)); 
		}
		}
	}
	rot.Work += tempW; // Diagnostics
}
