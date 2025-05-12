/** 
This page contains all functions associated with the fan. Note that this is the 3D extension off the well documented [fan page](http://basilisk.fr/sandbox/vheusinkveld/fan.c) 
*/

#include "fractions.h"	// Needed to compute fan fractions

struct sRotor rot;  		// Rotor details structure 
scalar fan[];			// Fan volume fraction

struct sRotor {
    bool fan;
    double start;		// Starting time
    double stop;		// Stopping time	
    double rampT;		// Time to start up rotor
    double P, Prho;		// Power, powerdensity 
    double R, W, A, V;		// Diameter, Thickness, Area ,Volume
    double diaVol;		// Diagnosed rotor volume 
    double x0, y0, z0;		// Origin of rotor
    double xt, yt, zt;      	// Translation angles
    double theta, phi;		// Polar and Azimuthal angle 
    double thetat, phit;    	// Translation angles 
    double Work;            	// Work done by rotor
    double cu;              	// Characteristic velocity
    bool rotate;            	// Rotation yes or no
    coord nf, nr;		// Normal vector fan, rotation 
};

/** Functions */
void init_rotor();
void rotor_update();
void rotor_coord();
void rotor_forcing();

/** Function returning the sRotor structure, includign default properties */
void init_rotor() {
    rot.Work = 0.;
    if(!rot.start)
    	rot.start = 0.;   
    if(!rot.stop)
    	rot.stop = 1E10;
    if(!rot.rampT)
    	rot.rampT = 30.;
    if(!rot.R)
	rot.R = 2.5;     
    if(!rot.W)
	rot.W = 0.3;    
    if(!rot.Prho)                  
     	rot.Prho = 2000.;		
    if(!rot.x0)
    	rot.x0 = L0/2.;
    if(!rot.y0)
	rot.y0 = 10.5;
    if(!rot.z0){
        #if dimension == 2
            rot.z0 = 0.;
        #elif dimension == 3
            rot.z0 = L0/2.;
        #endif
    }
    if(!rot.theta)
    	rot.theta = 97*M_PI/180.;		// Polar angle
    if(!rot.phi)
    	rot.phi = 0.*M_PI/180.;		// Azimuthal angle 

    if(rot.rotate) {
	if(!rot.phit){
	    rot.phit = 2*M_PI/240;
 	}
	if(!rot.xt){
	    rot.xt = 0;
 	}
       	if(!rot.yt){
	    rot.yt = 0;
 	}
	if(!rot.zt){
	    rot.zt = 0; 	
 	}
	if(!rot.thetat){
	    rot.thetat = 0.;
 	}

	} else {
        rot.xt = 0;
        rot.yt = 0;
       	rot.zt = 0;
       	rot.thetat = 0.;
       	rot.phit = 0.;
    }
    
    rotor_update();
}

/** Forcing by the rotor */
event forcing(i++) {
    if(rot.fan && t>rot.start && t<rot.stop) {
	rotor_coord();
	rotor_forcing();
    }
}

/** Rotate the rotor */
event rotate(i++) {
    if(t>rot.start && t<rot.stop){
        // Change center  
        rot.x0 += rot.xt;
        rot.y0 += rot.yt;
        rot.z0 += rot.zt;

        // Change angles 
        rot.theta += dt*rot.thetat;
        rot.phi += dt*rot.phit;

        rotor_update();
    }
}

/** Updating relevant rotor variables */
void rotor_update() {
    	rot.nf.x = sin(rot.theta)*cos(rot.phi);
	rot.nf.z = sin(rot.theta)*sin(rot.phi);
	rot.nf.y = cos(rot.theta);

	rot.nr.x = sin(rot.theta)*cos(rot.phi);
   	rot.nr.z = sin(rot.theta)*sin(rot.phi);
    	rot.nr.y = cos(rot.theta);

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
    	rot.cu = pow(3*rot.Prho*rot.W*crho, 1./3.);
}

/** Function returning the volume fractions of a fan object */
void rotor_coord() {
    scalar sph[], plnu[], plnd[];
    fraction(sph, -sq((x - rot.x0)) - sq((y - rot.y0)) - sq((z - rot.z0)) + sq(rot.R));
    fraction(plnu,  rot.nr.x*(x - rot.x0) + rot.nr.y*(y - rot.y0) + rot.nr.z*(z - rot.z0) + rot.W/2.);
    fraction(plnd, -rot.nr.x*(x - rot.x0) - rot.nr.y*(y - rot.y0) - rot.nr.z*(z - rot.z0) + rot.W/2.);	

    foreach () {
        fan[] = sph[]*plnu[]*plnd[];
    }
    boundary({fan});
}

/** Function returning new velocities based on a rotor forcing.
This is a function of powerdensity, width, direction, ramp-time, and diagnosed volume */

void rotor_forcing(){
    double tempW = 0.;
    double tempVol = 0.;
    double w, wsgn, damp, usgn, utemp, corrP;

    foreach(reduction(+:tempVol)){
         tempVol += dv()*fan[];
    }
    rot.diaVol = 1*tempVol;    
    Point point = locate(rot.x0, rot.y0, rot.z0);
    if(fan[] == 0.){
	foreach_dimension(){    	
             //Work in respective direction 
	     wsgn = sign(rot.nf.x*u.x[]) + (sign(rot.nf.x*u.x[]) == 0)*sign(rot.nf.x);
	     damp = rot.rampT + rot.start > t ? (t-rot.start)/rot.rampT : 1.;
	     w = wsgn*damp*sq(rot.nf.x)*(2.)*rot.P/pow(Delta,3)*dt;
	     tempW += fabs(w)/2*pow(Delta,3);
	
	     // New kinetic energy
	      utemp = sq(u.x[]) + w;
	      usgn =   1.*(u.x[] >= 0)*(utemp > 0) +
	              -1.*(u.x[] >= 0)*(utemp < 0) +
	     	       1.*(u.x[] <  0)*(utemp < 0) +
	     	      -1.*(u.x[] <  0)*(utemp > 0); 

	      //u.x[] = usgn*sqrt(fabs(utemp));
              u.x[] = usgn*min(sqrt(fabs(utemp)), sq(damp)*1.5*rot.cu); // Limiting the maximum speed

         }
     } else {

     foreach(reduction(+:tempW)) {	
	        if(fan[] > 0.) {
                foreach_dimension() {
	        // Work in respective direction 
		wsgn = sign(rot.nf.x*u.x[]) + (sign(rot.nf.x*u.x[]) == 0)*sign(rot.nf.x);
  	        damp = rot.rampT + rot.start > t ? (t-rot.start)/rot.rampT : 1.;
		corrP = rot.diaVol > 0. ? rot.V/rot.diaVol : 1.;
		w = wsgn*fan[]*damp*sq(rot.nf.x)*(2./rho[])*(corrP*rot.P/rot.V)*dt;
		tempW += dv()*fabs(w)/2;

		// New kinetic energy
		utemp = sq(u.x[]) + w;
		usgn = 	  1.*(u.x[] >= 0)*(utemp > 0) +
		    	 -1.*(u.x[] >= 0)*(utemp < 0) +
			  1.*(u.x[] <  0)*(utemp < 0) +
			 -1.*(u.x[] <  0)*(utemp > 0); 

		//u.x[] = usgn*sqrt(fabs(utemp));
		u.x[] = usgn*min(sqrt(fabs(utemp)), sq(damp)*1.5*rot.cu); // Limiting the maximum speed
	     }
	}
    }    
    }
    rot.Work += tempW;
}

// engine heat 10 percent efficiency, injected below the fan at 1 m above the surface.
event engine(i++) {
    double engine_heat = rot.P/0.1;
    double damp = rot.rampT + rot.start > t ? (t-rot.start)/rot.rampT : 1.;
    Point point = locate(rot.x0, 1., rot.z0);
    b[] += damp*engine_heat/(CP*pow(Delta,3))*gCONST/TREF*dt;
}
