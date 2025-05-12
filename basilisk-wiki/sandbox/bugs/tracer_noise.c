

/**
# Possible bug with tracer + embedded boundaries + two phases (miscible)*/

/**
There is a strange behaviour where "noise" appears on the tracer field when a solid is added to this scenario. */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};

/**
This following part was taken from two-phase-generic.h, and only the portion related to the viscosity computation on faces is included. If the entire file instead, convergence issues appear immediately (when embedded boundaries are present). This could be a separate issue or not (seems related to alpha field change from what I've seen).  */

#if 1
double mu1 = 0., mu2 = 0.;

event defaults (i = 0)
{
	if (mu1 || mu2)
		mu = new face vector;
}

#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
#endif

#define sf f
event properties (i++)
{
	foreach_face() {
		double ff = (sf[] + sf[-1])/2.;
		if (mu1 || mu2) {
			face vector muv = mu;
			muv.x[] = fm.x[]*mu(ff);
		}
	}
}
#endif

//#include "two-phase-generic.h"


/**
The setup for the scenario is below: */

int maxlevel = 5;

//Parameters
double U1 = 0.001;
double H1 = 0.000249999;
double tend =  1.5;

double f_error = 1e-2;
double u_error;


u.n[left]  = dirichlet(U1*cs[]);
u.t[left]  = dirichlet(0.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

f[left] = 1.;

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

/**
A boolean variable embSolid is used to execute the case with and without the solid:
*/

bool embSolid = true;
int main()
{
	L0 = 5.*H1;
	origin (0.,-L0/2.);
	N = 64;
	DT = 0.01;
	u_error = 0.05*U1;

	mu1 = 0.0966;
	mu2 = 0.0023;
	
	run();
	embSolid = false;
	run();
}

event end (t=tend)
{}

event init (t = 0)
{

/**
If the solid is added, the "noise" on the tracer is present (also the statistic mgu.i for the viscosity grows considerably). Otherwise the effect is not present: */
	if (embSolid){
		solid (cs, fs,  intersection(-y+H1*2. ,+y+H1*2.));
		fractions_cleanup (cs, fs, 0.0001);
	}
}

/**
The adaptive mesh refinement used is below. Completely removing adaptive mesh refinement also makes the effect disappear: */

event adapt (i++) {
	if(embSolid)
		adapt_wavelet ({cs,u,f}, (double[]){1e-3,u_error,u_error,f_error}, maxlevel);
	else
		adapt_wavelet ({u,f}, (double[]){u_error,u_error,f_error}, maxlevel);
}

/**
Output generation: */

event output (t = 1.2){
	if(embSolid)
		output_ppm(f, file = "tracer_noise_1.png", n = 400, map=cool_warm, spread=3);
	else
		output_ppm(f, file = "tracer_noise_2.png", n = 400, map=cool_warm, spread=3);
}

/**
#Results obtained in both cases
The results obtained for each setup are shown below:

![This image shows how the tracer is affected and a "noise" appears when the solid is present.](tracer_noise/tracer_noise_1.png)

![This image shows the result obtained when the solid is not present and the tracer diffuses freely.](tracer_noise/tracer_noise_2.png)

*/