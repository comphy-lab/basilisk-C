/**
# Wave generation by fluidised flow. Three-phase Navier-Stokes with option for diffusion between two phases. This file implements a setup for a fixed input condition, 10 cm from the impact zone, where $h_f$ (impact height), $u_f$ (impact velocity) and $\theta$ can be defined by the user (UF, H0 AND THETA0). Please see [example script](http://basilisk.fr/sandbox/pdc_final/script).

This is the Navier-Stokes VOF case of the [Bougouin, 2020](#references) fluidized granular flow tsunami generation experiment.
*/
#include <sys/stat.h>
#include "grid/quadtree.h"
#include "navier-stokes/centered.h"

scalar f3[];
double rho3 = 1400; //density of dense salt water
double mu3  = 0.1; //density of salt water at 20degrees

#define THETA0 0.262
#define BHO 0.02
#define H0 0.02
#define UF 2.2


#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(clamp(f3[],0.,1.)*(rho3 - rho1) + rho1 -rho2) + rho2)
#endif
#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(clamp(f3[],0.,1.)*(mu3 - mu1) + mu1 - mu2) + mu2) 
#endif


#define MAXLEVEL 12

/** Definition of the robin boundary, see [Antoon's robin.c](http://basilisk.fr/sandbox/Antoonvh/robin.c) for details. */
#define robin(a,b,c) ((dirichlet ((c)*Delta/(2*(b) + (a)*Delta))) + ((neumann (0))* ((2*(b) - (a)*Delta)/(2*(b) + (a)*Delta) + 1.)))

#include "two-phase.h"
#include "../myconserving.h" 
#include "utils.h"

double ue = 1e-3;
FILE * fp;

/** Define parameters */
#define io (0.01) // Width of initial jet
#define sl (0.11) // Length of exposed slope
#define add_length (0.265/sin(THETA0)) // 0.265 IS THE WATER DEPTH (M)
#define MUR 54
#define imagerate 0.01
#define endlim 0 // UPDATE
#define TRIANGLE (tan(THETA0)*(x-(io+sl+add_length)))

int main() {
	  L0 = 12.;
	    origin (0, 0);
	    rho1 = 998; //water
	    rho2 = 1.27; //air
	
	    mu1 = 1.002e-3; //water
            mu2 = 1.002e-3/MUR; //air	
  
	    struct stat st = {0};
	    if (stat("./profiles", &st) == -1 && stat("./images", &st) == -1)  {
		    mkdir("./profiles", 0755);
		    mkdir("./images",0755);
		    
	    }
	    
	    init_grid (1 << (MAXLEVEL-1));
	    fp=fopen("tracers_cons_final.dat", "w");
 	    const face vector av[] = {9.81*sin(THETA0), -9.81*cos(THETA0), 0};
	    system("mkdir el");
	    a = av;

	    CFL = 0.3;
	    DT=0.005;
	    run(); 
}


/** Set boundary conditions */

u.t[bottom] = robin(1.,BHO,0.); //robin BC - 0.015 should be the height of the 
u.n[left] = dirichlet((y < H0 && t < 0.5) ? UF : 0 );
f3[left] = dirichlet(y < H0 && t < 0.5);


/** Initilaise domain */

event init (t = 0)
{
	 
	  refine(y > 3 && level < 7);
	  foreach() {
		  if (y < tan(THETA0)*(x-(sl+io))) { //with water
			  f[] = 1.;
			  f3[] = 0.;
		  }
		  else if (y <= H0 && x < io) {
			  f3[] = 1.;
			  f[] = 1.;
		  }
		  else {
			  f[] = 0.;
		          f3[] = 0.;
		  }
		  foreach_dimension()
			  u.x[] = 0.;
		}
	  boundary({f, f3, u});

}

/** Slope implementation */
scalar tri[];
event make_slope(i++) {
	fraction (tri, y <= TRIANGLE);
	foreach()
		foreach_dimension()
			u.x[] -= u.x[]*tri[];
  /*
	face vector tri_face[];
	face_fraction (tri, tri_face);
	boundary((scalar *){tri_face});
	foreach_face()
		uf.x[] -= tri_face.x[]*uf.x[];
                */
	boundary ((scalar *){u,uf});
}

/** Log files */
#include "view.h"
event logfile (i++){
	  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}


/** Output files */
int time_step = 0.;
double shift_control = 1.;
event output_files  (t += imagerate; t < endlim) { 
                        // Images (to postprocess in python)
                        char one[80], two[80], three[80], four[80];
			foreach() {
		     		if (tri[] > 0.5)
					tri[] = -1;
			}		
			sprintf (two, "images/density-%d.png",time_step );  
			output_ppm(rho, file = two, min = 1, max = 1400, n = 2048, linear = true, box = {{0,0},{3,1}}, mask = tri);
			
                        // Velocity profiles
			sprintf (three, "profiles/vprof-%d",time_step);
			FILE * fp1 = fopen (three, "w");
			//double step = (L0)/(1 << (lmax));
			for (double y = 0.; y <= 0.06; y += 8/(pow(2,MAXLEVEL))){
				fprintf (fp1, "%g %g %g %g\n", y,interpolate (u.x, 0.02, y),interpolate (u.y, 0.02, y),interpolate(f, 0.02,y));} //Outputs y, u.x, u.y, f
			fclose (fp1);

			time_step += 1.;
}


/** Output VOF facets */
int add = 0;
event interface (t+= 0.004) {
		char fname[99];
		sprintf(fname, "el/water_elevation%d.txt", add);
		FILE * fp = fopen(fname, "w");
		output_facets(f,fp);
		fclose(fp);
		add += 1.;

}

/** Energy calculations */
int dissipation_rate (double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  double rateFlow = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir) reduction (+:rateFlow)) {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			     sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			     sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;

    rateWater += mu1/rho1*(1.-f3[])*f[]*sqterm; //water (when f3 is one, this will be zero)
    rateFlow += mu3/rho3*(f3[])*f[]*sqterm; //flow
    rateAir   += mu2/rho2*(1. - f[])*sqterm; //air
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  rates[2] = rateFlow;
  return 0; }



event graphs (t += 0.01) {
  static FILE * fpwater = fopen("budgetWater.dat", "w");
  static FILE * fpair = fopen("budgetAir.dat", "w");
  static FILE * fflow = fopen("budgetFlow.dat", "w");

  double keWater = 0., gpeWater = 0.;
  double keAir = 0., gpeAir = 0.;
  double keFlow = 0., gpeFlow = 0.;
  foreach(reduction(+:keWater) reduction(+:gpeWater) 
	  reduction(+:keAir) reduction(+:gpeAir)
	  reduction(+:keFlow) reduction(+:gpeFlow)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    keWater += rho1*(1.0-f3[])*norm2*f[]*dv();
    keAir += rho2*norm2*(1.0-f[])*dv();
    keFlow += rho3*f3[]*norm2*f[]*dv();

    gpeWater += rho1*f[]*(1.0-f3[])*dv()*9.81*(-sin(THETA0)*x+cos(THETA0)*y);
    gpeAir += rho2*(1.0-f[])*dv()*9.81*(-sin(THETA0)*x+cos(THETA0)*y);
    gpeFlow += rho3*f[]*(f3[])*dv()*9.81*(-sin(THETA0)*x+cos(THETA0)*y);

  }
  double rates[3];
  dissipation_rate(rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
  double dissFlow   = rates[2];

    if (i == 0) {
    fprintf (fpwater, "t ke gpe dissipation\n");
    fprintf (fpair, "t ke gpe dissipation\n");
    fprintf (fflow, "t ke gpe dissipation\n");
    }
  fprintf (fpwater, "%g %g %g %g\n",
	   t, keWater/2., gpeWater, dissWater);
  fprintf (fpair, "%g %g %g %g\n",
	   t, keAir/2., gpeAir, dissAir);
    fprintf (fflow, "%g %g %g %g\n",
	   t, keFlow/2., gpeFlow, dissFlow);

}



/** Adaption */
event adapt (i++)
{
	scalar omega[];
	vorticity(u,omega);
	 //adapt_wavelet((scalar *){cs,u,f},(double[]){ue,ue,ue,0.1},maxlevel,(maxlevel-4));
	adapt_wavelet((scalar *){omega,f},(double[]){ue,ue,0.01},MAXLEVEL-1);
				        
}
