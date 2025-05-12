/**
# Wave generation by fluidised flow. Three-phase Navier-Stokes with option for diffusion between two phases.

This is the Navier-Stokes VOF case of the [Bougouin, 2020](#references) fluidized granular flow tsunami generation experiment.

*/

//#include "grid/octree.h"
#include "navier-stokes/centered.h"

/** Initialise density tracer f3 (representative of third phase)*/
scalar f3[];
double rho3 = 1400; //density of dense salt water
double mu3 = 0.1; //density of salt water at 20degrees

/**
  We “overload” the default by defining the rho() and mu() macros before including the code for two phases. */

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(clamp(f3[],0.,1.)*(rho3 - rho1) + rho1 -rho2) + rho2) 
#endif
#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(clamp(f3[],0.,1.)*(mu3 - mu1) + mu1 - mu2) + mu2) //viscoisty is a function of the volumne fraction
#endif

/** Definition of the robin boundary, see [Antoon's robin.c](http://basilisk.fr/sandbox/Antoonvh/robin.c) for details. */
#define robin(a,b,c) ((dirichlet ((c)*Delta/(2*(b) + (a)*Delta))) + ((neumann (0))* ((2*(b) - (a)*Delta)/(2*(b) + (a)*Delta) + 1.)))


/**Optional to include diffusion */
#def IMPLICIT_DIFFUSION 0

#if IMPLICIT_DIFFUSION
#include "sandbox/qmagdelaine/phase_change/advection_Q.h" 
#endif

#include "two-phase.h"
#include "myconserving.h" //updated conserving for new three-phase formulation

#if IMPLICIT_DIFFUSION
#include "sandbox/qmagdelaine/phase_change/mixtures.h"
//#include "sandbox/qmagdelaine/phase_change/elementary_body.h" //Optional
#endif
//#include "tension.h" //Optional: if including surface tension
#include "utils.h"

double ue = 1e-3;

#define io (0.338)
#define MUR 54
#define MAXLEVEL 12

/** Initialise domain, 6 units long. */
int main() {
	  L0 = 6.;
	  origin (0, 0);
	  rho1 = 998; //water
	  rho2 = 1.27; //air
	  mu1 = 1.002e-3; //water
          mu2 = 1.002e-3/MUR; //air	    
	  //f.sigma = 73e-3;
	  
	  init_grid (1 << MAXLEVEL);
	  //fp=fopen("tracers_cons_final.dat", "w");
 	  const face vector av[] = {9.81*sin(0.268), -9.81*cos(0.268), 0};
	  a = av;
	  CFL = 0.4;
	  run(); 
}

/**Implementation of robin boundary condition, robin(a,b,c) where b is the slip length. */
u.t[bottom] = robin(1.,0.04,0.); 
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);


/** Optional if including diffusion */
#if IMPLICIT_DIFFUSION
attribute {
	  double D;
}
#endif

/** Initialise the geometry. */
event init (t = 0)
{
	  foreach() {
		  if (y < tan(0.268)*(x-(1+io))) { //with water
			  f[] = 1.;
			  f3[] = 0.;
		  }
		  else if (x < ((y - tan(1.309)*(io))/(-tan(1.309))) && y < 0.225 && z <= 0.2 ) { //accounts for the extra volume
			  f[] = 1.;
			  f3[] = 1.;
		  }
			  else {
			  f[] = 0.;
		          f3[] = 0.;
		  }
		  foreach_dimension()
			  u.x[] = 0.;
		}
	boundary ({f,f3,u});

}

//* Output statistics to a logfile. */
#include "view.h"
event logfile (i++){
	  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}


/**
If accounting for diffusion... fixme: diffusion into bottom boundary!*/
#if IMPLICIT_DIFFUSION
#define D_L 0.0001
mgstats mgd1;
event tracer_diffusion (i++) {
	  foreach()
                f[] = clamp(f[], 0., 1.);
	  boundary ({f});

	  foreach()
              f3[] = (f[] > F_ERR ? f3[]/f[] : 0.);
	  boundary({f3});
          f3.D = D_L;
          f3.inverse = false;
          no_flux_diffusion (f3, f, dt);
	  foreach()
		f3[] *= f[];
	  boundary({f3});
			

}
#endif

/** Crude method for boundary implementation but works well for no diffusion case:*/
scalar tri[];
#define TRIANGLE (tan(0.268)*(x-(1+io+1.02)))
event makeslope (i++) {
                fraction (tri, y <= TRIANGLE);
                foreach()
			foreach_dimension() 
				u.x[] -= u.x[]*tri[];
		face vector tri_face[];
 		face_fraction (tri, tri_face);
 		boundary ((scalar *){tri_face});
 		foreach_face()
                        uf.x[] -= tri_face.x[]*uf.x[];
  		boundary ((scalar *){u, uf});  

}

//3D boundary
scalar edge[];
event makeedge (i++) {
	fraction(edge, z > 0.2);
	foreach()
		foreach_dimension()
			u.x[] -= u.x[]*edge[];
	boundary((scalar *){u, uf});
}



/**Output images*/

int time_step = 0.; //don't need this
event output_files  (t += 0.01; t < 4.5) { 
			
			char one[80], two[80], three[80], four[80];
  
			foreach() {
				   if (tri[] > 0.5)
					       tri[] = -1;
                        }
		     	   
			//First two outputs show JPG snapshots for FOV 
			sprintf (one, "tracer-%d.png",time_step );
			output_ppm(f3, file = one, map = cool_warm, n = 2048, linear = true, box = {{0,0},{4.5,1}}, mask = tri);
			sprintf (two, "density-%d.png",time_step );  
			output_ppm(rho, file = two, min = 1, max = 1400, n = 2048, linear = true, box = {{0,0},{4.5,1}}, mask = tri);
			
  
  
                        //Bview outputs
			view(fov = 20, psi = 0.268, tx = -0.55, ty = -0.05);
			box();
			char str[80];
			sprintf(str, "t = %g", t);
			draw_string(str, pos = 4, lw =3);
			draw_vof("f", filled = -1, fc = {1,1,1});
						squares("f3");
			save("movie.mp4");
			
			//We also output a vertical cross section of the velocity compoentns 
			sprintf (three, "vprof-%d",time_step);
			FILE * fp1 = fopen (three, "w");
			for (double y = 0.; y <= 0.05; y += 8/(pow(2,MAXLEVEL))){
				fprintf (fp1, "%g %g %g %g\n", y,interpolate (u.x, 1.238, y, 0.1),interpolate (u.y, 1.238, y, 0.1),interpolate(f, 1.238,y, 0.1));} //Outputs y, u.x, u.y, f
			fclose (fp1);
			time_step += 1.;
}

/**
## Outputs

We are interested in the viscous dissipation rate in both water, granular fluid and air. */

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

/**
We log the evolution of the kinetic and potential energies and
dissipation rate as functions of the non-dimensional time. */

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

    gpeWater += rho1*f[]*(1.0-f3[])*dv()*9.81*(-sin(0.268)*x+cos(0.268)*y);
    gpeAir += rho2*(1.0-f[])*dv()*9.81*(-sin(0.268)*x+cos(0.268)*y);
    gpeFlow += rho3*f[]*(f3[])*dv()*9.81*(-sin(0.268)*x+cos(0.268)*y);

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


/**
We output the water elevation profile, which is used to extract wave gauges. */
int add = 0;
event interface (t+= 0.01) {
		char fname[99];
		sprintf(fname, "water_elevation%d.txt", add);
		FILE * fp = fopen(fname, "w");
		output_facets(f,fp);
		fclose(fp);
		add += 1.;

}

/**Option not to include - vtu files not necessary.  */ 

int counter = 0;
#include "output_vtu_foreach.h"
event vtk_file (t += 0.05) {
			FILE * fp;
			char name[80];
			sprintf(name, "test_0.01_%d.vtu", counter);
			fp = fopen(name,"w");
			output_vtu_bin_foreach((scalar *) {rho,f}, (vector *) {u},N,fp,false);
			fclose(fp);
			counter += 1.;
}




/**Adapt with respect to vorticity. */ 

event adapt (i++)
{
	scalar omega[];
	vorticity(u,omega);
	adapt_wavelet((scalar *){omega,f},(double[]){ue,ue,0.01},MAXLEVEL);
				        
}


/**
## References
A. Bougouin,  R. Paris and O. Roche. Impact of Fluidized Granular Flows into Water: Implications for Tsunamis Generated by Pyroclastic Flows. 
*/