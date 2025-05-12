/** 
This page contains al physics accosciated functions of the Idealized case. Note that the buoyancy implementation is similar to [here](http://basilisk.fr/sandbox/Antoonvh/vhm16.c) and the lumped soil model is documented [here](http://basilisk.fr/sandbox/Antoonvh/lumpedsoil.c). Note that for the imported [SGS.h](http://basilisk.fr/sandbox/Antoonvh/SGS.h) page, the diffusion event is turned off since this is implemented here. 

What is done: 
Profiles are set. Boundary conditions and profiles are initialized. The vreman subgrid closure is used. Gravity acceleration is added via buoyancy. A lumped soil model is used, tendency towards a soil temperature. 
*/

#if dimension == 3
	#include "SGS.h"
#endif 

#define CP 1005.	// C_p for air 
#define gCONST 9.81	// Gravitational constant
#define TREF 273.	// Kelvin
#define INVERSION .3 	// Kelvin per meter
#define karman 0.4      // von Karman constant 

#define roughY0u 0.1    // roughness wind length 
#define roughY0h 0.1     // roughness heat length

#define WIND(s) 0. //(-max(0.2*log(s/roughY0u),0.))   // log Wind profile 
#define QFLX 0. 	// 0 (0.001 = 20wm2)
#define BSURF (1.5*b[] - 0.5*b[0, 1])
#define GFLX (-Lambda*(BSURF - bd))
double Lambda = 0.005, bd = 0.;   // Grass coupling
#define STRAT(s) gCONST/TREF*s*INVERSION

scalar b[];
scalar * tracers = {b};

double crho = 1.;	
face vector av[]; 
struct sCase def;

struct sCase {
	double wind;
	double wphi;
};	

void init_physics(){
 	def.wind = 1.;
        def.wphi = 0.;

	b.nodump = false; 
        u.n[bottom] = dirichlet(0.);
        u.t[bottom] = dirichlet(0.);
        u.n[top] = dirichlet(0.);
        u.t[top] = dirichlet(WIND(y));

        periodic (left);
	
	b[bottom] = BSURF;
	b[top] = dirichlet(STRAT(y));
	
	#if dimension == 3
            u.r[top] = dirichlet(WIND(y));
	    u.r[bottom] = dirichlet(0.); 
            u.t[top] = neumann(0.); 
		
	    Evis[bottom] = dirichlet(0.); // Flux is explicitly calculated
            Evis[top] = dirichlet(0.);

	    periodic(front);
	#endif  
	foreach() {
		b[] = STRAT(y);
		if(fabs(def.wind) > 0.) {
       	            u.x[] = WIND(y);
		}
	}
}
double lut[20];
double lut2[20];
event init (t = 0) {
    for (int m = 0; m <= 19; m++) {
	double d  = (L0/((double)(1 << m)))/roughY0u;
	double d2 = (L0/((double)(1 << m)))/roughY0h;

        if (m==0) {
	    lut[0] = 0;
	    lut2[0] = 0;
	} else {
            lut[m] = sq(karman/(log(d) - 1.));
	    lut2[m] = (log(4.*d2) - 1.)/(log(d2) - 1.);
        }
    }
}

/* Gravity forcing */
event acceleration(i++){
	foreach_face(y){
		av.y[] = (b[] + b[0,-1])/2.;
	}
}

event inflow(i++){
    double sides = 50;
    double relaxtime = dt/50;
    foreach(){
	if((x < sides || x > L0-sides) ||
	   (z < sides || z > L0-sides) ||
	   (y > L0-2*sides )) {
	    double a = (x < sides) ? x : fabs(x-L0);
	    a = 1.; 
	    u.x[] = u.x[] + a*(WIND(y)-u.x[])*relaxtime;
 	    b[] = b[] + a*(STRAT(y) - b[])*relaxtime;
	    u.y[] = u.y[] - a*u.y[]*relaxtime;
            u.z[] = u.z[] - a*u.z[]*relaxtime;
	}
    }
}

mgstats mgb;
/* Diffusion */
event tracer_diffusion(i++){
    scalar r[];
    foreach() {
        r[] = 0;
        if (y < Delta)
            r[] = (QFLX + GFLX)/sq(Delta); // div needed as normalization 
    }
    
    double flx = 0, bt = 0;
    double fctr = CP*TREF/gCONST;
    foreach_boundary(bottom reduction(+:flx) reduction(+:bt)) {
        flx = flx + (QFLX + GFLX) * sq(Delta);
         bt = bt + BSURF * sq(Delta);
    }
    bt = bt/sq(L0);
    flx = flx/sq(L0);
    fprintf(stderr, "soil=%g %g %g %d\n", t, fctr*flx, fctr*bt/CP, i);  
    
    mgb = diffusion(b, dt, mu, r = r);
}
