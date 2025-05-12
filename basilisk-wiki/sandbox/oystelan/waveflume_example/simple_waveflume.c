
/** 
# Simple flume with wavemaker
This little example which shows how to implement a simple waveflume in Basilisk, where a wave is propagated into the domain from the left boundary. To keep it simple, linear wave theory is used. The example may however easily be replace by more advanced theories by replacing the functions velox, veloy and wave_elev.
The wave elevation at inflow is measured using the waveprobe functionality added with "waveprobes.h".

Oystein Lande 2018

We start by including the building blocks:
*/
#include <sys/stat.h>
#include "grid/quadtree.h"
#include "adapt_wavelet_leave_interface.h"
#include "utils.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "waveprobes.h"
#include "output_vtu_foreach.h"
#include "reduced.h"
#include "tension.h"
#include "view.h"



/** For simplicity, all input parameters defining the simulation is defined in the struct below */

struct SIM_properties {
	double tmax; // simulation length
	double dtmax;
	double size[4]; // wave flume dimensions, [x_min, x_max, y_min, y_max]
  	double waterlevel; // elevation of stillwater level, relative to domain definition
  	double gravity; // m/s^2 (default 9.81)
  	double rho1; // density of water
  	double rho2; // density of air
    int refine_init; // refinment level at t = 0
	int adapt_refine_max; // maximum level of refinement
	int adapt_refine_min; // minimum level of refinement
	int adapt_utol;
	double A; // wave linear amplitude
	double w; // wave frequency (rad/s)
	double k; // wave number
	double waterdepth; // waterdepth
};

struct SIM_properties simdata = {
	.tmax = 20.0,
	.dtmax = 0.5,
	.size = {-200.,200.,-90.0,30.0},
	.waterlevel = 0.0,
	.gravity =9.81,
	.rho1 = 1025,
	.rho2 = 1.225,
    .refine_init = 9,
	.adapt_refine_max = 9,
	.adapt_refine_min = 4,
	.adapt_utol = 1.0,
	.A = 0.25,
	.w = 1.0472,
	.k = 0.11179,
	.waterdepth = 90.
};

int vtucount = 0;

/**
## Functions
We define some useful functions which will be called further down.

Sets the domain size and origin according to the specified values in simdata
*/
void set_domain_size(){

	double lx = simdata.size[1]-simdata.size[0];
	double lz = simdata.size[3]-simdata.size[2];

	fprintf(stderr,"data:%.2f\n",simdata.size[1]);
	init_grid (1 << (simdata.refine_init));

	if (lx >= lz){
		size (lx);
		origin (simdata.size[0], simdata.size[2]); // move origin
		//mask (y > dp.size[5]  ? top : none);
	}
	else {
		size (lz);
		origin (simdata.size[0], simdata.size[2]); // move origin
		//mask (x > dp.size[1]  ? right : none);
	}
}


/** maskes away parts of the domain to make it the right size*/
void  mask_domain(){

	double lx = simdata.size[1]-simdata.size[0];
	double lz = simdata.size[3]-simdata.size[2];

	if (lx >= lz){
		mask (y > simdata.size[3]  ? top : none);
	}
	else {
		mask (x > simdata.size[1]  ? right : none);
	}
}

/** fills the basin and sets kinematics to 0 in entire basin*/
void set_kinematics(){
	// initialize basin with wave
	scalar phi[];	
	foreach_vertex() {
		phi[] = -y + simdata.waterlevel;
	}
    fractions (phi, f);
	fprintf(stderr,"Initializing basin velocities... rest...\n");
	foreach()
		foreach_dimension()
	        u.x[] =    0.0;
	boundary((scalar*){f,u});
}

/** Nice to have print function*/
void mg_print (mgstats mg)
{
    if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 exp (log (mg.resb/mg.resa)/mg.i));
}

/** horizontal particle velocity, linear wave theory for finite depth*/
double velox(double x,double y,double t){
		return clamp(t/4.,0.,1.)*simdata.A*simdata.w*(cosh(simdata.k*(y+simdata.waterdepth))/sinh(simdata.k*simdata.waterdepth))*cos(simdata.k*x - simdata.w*t);
}

/** Vertical particle velocity, linear wave theory for finite depth*/
double veloy(double x,double y,double t){
		return clamp(t/4.,0.,1.)*simdata.A*simdata.w*(sinh(simdata.k*(y+simdata.waterdepth))/sinh(simdata.k*simdata.waterdepth))*sin(simdata.k*x - simdata.w*t);
}

/** wave elevation, linear wave theory*/
double wave_elev(double x, double t){
		return clamp(t/4.,0.,1.)*simdata.A*cos(simdata.k*x - simdata.w*t);
}

/** Simple coarse function for calculating volume fraction */
double wave_Vfrac(double x, double y, double t, double del){

	double eta = wave_elev(x,t);
	//cout << wwelev << endl;
	if (eta < y - (del / 2.)) {
		return 0.0;
	}
	else if (eta > y + (del / 2.)) {
		return 1.0;
	}
	else {
		// Calculate volume fraction for the given cell with size del and position y
		return (eta - (y - (del / 2.))) / del;
	}
}

/** write unstructured vtu files */
void save_vtu ( int nf, int j)
{
  char name[80];
  FILE * fp ;
  nf > 0 ? sprintf(name, "RES_VTK/res_n%3.3d_%4.4d.vtu",pid(),j) : sprintf(name, "RES_VTK/res_%4.4d.vtu",j);
  fp = fopen(name, "w"); output_vtu_bin_foreach ((scalar *) {f,p}, (vector *) {u}, N, fp, false); fclose (fp);
}


/** 
## Setting boundaries */
// Left boundary (wave inflow)
u.n[left]  = f[]*dirichlet(velox(x,y,t)) + (1.-f[])*neumann(0);
u.t[left]  = f[]*dirichlet(veloy(x,y,t)) + (1.-f[])*neumann(0);
f[left]    = wave_Vfrac(x,y,t,Delta);

// Top boundary
u.n[top] = neumann(0);
p[top] = dirichlet(u.y[]*abs(u.y[])*rho(f[])*0.5); // this is a neat little trick to avoid escalating back-circulating flows in the top boundary.

// bottom boundary
p[bottom]  = neumann(-a.y[]*fm.y[]*rho(f[])); // Think this is already default, but just to be sure...

/**
## Main loop
*/
int main() {
    TOLERANCE = 1E-8;
    mkdir("./RES_VTK",0755); // make a directory to store the resulting VTU files

	//Set domain properties
	set_domain_size();

    dtmax = simdata.dtmax;
    rho1 = simdata.rho1, rho2 = simdata.rho2;
    mu1 = 0., mu2 = 0.;

    G.y = -simdata.gravity;
    Z.y = 1.;

    run();

}


/**
### initialize basin at rest
*/
event init (t = 0) {
	// Mask domain to fit specified size
	mask_domain();
	// Set flume surface at rest
	set_kinematics();
}

/**
# Events
*/

/** This is a useful function to dump some simulation data during runtime, but strictly not neccessary */
#if 1
event logfile (i++) {
    stats s = statsf (f);
    scalar l[];
    foreach()
    l[] = (sqrt(sq(u.x[]) + sq(u.y[])));

    norm n = normf (l);
    fprintf (stderr, "time: %g i: %d dt: %g statsF: %g %g %g", t, i, dt, s.sum, s.min, s.max - 1.);
    fprintf (stderr, ", Urms: %g Umax: %g, speed: %g cellcount: %ld\n", n.rms, n.max, perf.speed, grid->tn);

    mg_print (mgp);
    //mg_print (mgpf);
    mg_print (mgu);
    fflush (stderr);
}
#endif


/** 
## Waveprobes
This function uses the height function to calculate the surface elevation at a given point x (or x,y in 3D)
*/
#if 1
event waveprobe (t+=0.02;t<=simdata.tmax) { 
	heights (f, h);
	static FILE * fp0 = fopen("waveprobe0.dat", "w");
	double ycoords[2]  = {-10.,10.}; // define a vertical line of points 
	double yMax0 = wprobe(-199.9,ycoords,20);
	double yMax1 = wave_elev(-199.9,t);
	// update file
	fprintf(fp0, "%g %g %g\n",t, yMax0,yMax1);
	fflush(fp0);
}
#endif


/** dump a vtu file (which can be viewed directly in paraview), for a given timeinterval */
#if 1
event logfilevtu (t=0.0;t<=simdata.tmax;t+=0.04) {
  save_vtu(0,vtucount);
  vtucount += 1;
}
#endif


/**
## mesh adaptation
Using a sligtly modified version of adapt_wavelet, which restricts the mesh to be at max level around the fluid interface
*/
#if 1
event adapt(i++){
	adapt_wavelet_leave_interface((scalar *){u},{f},(double[]){simdata.adapt_utol,simdata.adapt_utol}, simdata.adapt_refine_max, simdata.adapt_refine_min,1);
}
#endif

/**
#Results
The simulation takes about a min on a resonable desktop computer. 


*/
