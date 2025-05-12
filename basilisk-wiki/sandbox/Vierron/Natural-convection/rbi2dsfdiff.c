/** ## Rayleigh Benard instability with stratified fluid 2D and Adaptative mesh */
#include "tag.h"
#include "view.h"
#include "convection_boussinesq_buoyancy-adapt.h"
#include "global_nusselt.h"
#include "profil5c.h"
#include "navier-stokes/perfs.h"
#include "curvature.h"

/* Level of refinement of the adaptive mesh size 2**MINLEVEL
   Here the adaptive mesh is disabled because grid/multigrid.h is used. */
#define MINLEVEL 6
#define MAXLEVEL 8
#define L0 4.

/* Initial height of the fluid layer. */
double A = 0.2; 

/* Max time step. */
double EndTime= 300.;
double ytop=0.5;
double ybot=-0.5;

/* The domain depends on the number of processors affected in order to be able to run the code in parallel.
   We have an aspect ratio of one processor at y and 4 at x. Gives it a DT variable 
   which is the maximum time step to help run the code. 
   TOLERANCE corresponds to the minimum value to be reached for residues (poisson.h). */
int main() {
	size (L0); 
	origin (0., -L0/2.);
	//DT = 0.1;
	//TOLERANCE = 1e-6;
	init_grid(1<<MAXLEVEL);
	Ra = 1e3; Pr = 1000.; B = 1.1;
	run();
}

/** 
## Boundaries Conditions 
*/

T[left] = neumann(0.);
T[right] = neumann(0.);

T[embed] = y > 0.2 ? dirichlet(-0.5) : dirichlet(0.5);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);


/** ## Initialization of the fluid layer + temperature field + velocity field. */
event init (t=0) {
	vertex scalar phi[];
        foreach_vertex() {
           phi[] = intersection (0.5 - y, 0.5 + y);
        }
        boundary ({phi});
        fractions (phi, cs, fs);

	fraction (f, - (y - Y0 - A));
  	foreach(){
   		T[] = - y;
    		foreach_dimension()
      			u.x[] = 0.;
 	}
	boundary ({T,u,f});
}

/** ## Adaptative mesh */
event adapt (i++) {
	adapt_wavelet ({T,f,u}, (double[]){1e-2,1e-3,1e-2,1e-2}, MAXLEVEL);
}

/** ## Output logfile
*/
/* scalar div corresponds to the field of divergence
   stats s0 allows to evaluate statistical quantities on this field
   s0.sum/s0.volume is the average value of the divergence
   s0.max is the maximum value of this field
   It is important to check these quantities because the method chosen 
   by Basilisk roughly forces the null divergence. */
event logfile (t += 1.0; t <= EndTime) {
	scalar div[];
	foreach() {
		div[] = 0.;
		foreach_dimension()
			div[] += u.y[1] - u.y[];
		div[] /= Delta;
	}
  /* statsf(div) function returns the minimum, maximum, volume sum, standard deviation and volume for field div.
     The mgT is the statistics on the Poisson.h solver for the diffusion step which is in the convection_boussinesq.h */
	stats s0 = statsf (div);
	static FILE * fs = fopen("stat","w");
	fprintf (fs, "%f %.9g %.9g %d %d\n",
		   t, s0.sum/s0.volume, s0.max, mgT.i, mgT.nrelax);
          fflush (fs);

/** ## Outputs
*/
  /* Creation of video file mp4 of the temperature field. */
  	output_ppm (T, file="temperature.mp4", n = 1024, box = {{0.,ybot},{L0, ytop}});
  	scalar l[];
	foreach(){
		l[] = level;
	}
	output_ppm (l, file="grid.mp4", n=512, min=MINLEVEL, max=MAXLEVEL, box = {{0.,ybot},{L0,ytop}});

  /* Variable initialization nusselt + nusselt calculation (global_nusselt.h). */
#if 1
	double nu_vol=0. , nu_t=0. ,nu_b=0.;
	nu_vol = nusselt_vol(T,u); 
	nu_t = nusselt_top(T); 
	nu_b = nusselt_bot(T);

  /* Creation of the data file in order to store my variables at each time step (t += 1.0; t <= EndTime). */
  	static FILE * fy = fopen ("data", "w");
	if (t<1.){
	    fprintf (fy, "#[1]Ra [2]Pr [3]N [4]nutop [5]nubot [6]nuvol [7]umin [8]umax [9]vmin [10]vmax [11]t\n");
 	 }
  /* statsf(velx and vely) returns the minimum, maximum, volume sum, standard deviation and volume for field velx and vely.*/
  	stats velx = statsf (u.x), vely = statsf (u.y);
  	fprintf (fy, "%1.1e %.1f %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g %f \n",
	   Ra,Pr,N,nu_t,nu_b,nu_vol, velx.min, velx.max, vely.min, vely.max, t);
  
  /* Creation of video file mp4 of the stratification. */
  	output_ppm (f, file="strat.mp4", n = 1024, box = {{0.,ybot},{L0, ytop}}, map = gray, min = 0., max= 1.);

  /* Creation of video file mp4 of the stratification + the temperature field. */
	clear();
	view(tx=-0.25);
	box();
	draw_vof("f");
	squares ("T", spread=1., linear = true);
	save("temp+strat.mp4");
}

/* Calculation of average temperature profile (Antoonvh Sandbox) + numerical data temperature field
   + saving of the simulation data set at t = EndTime (dump). */
event pro (t=EndTime) {
	profile({T,u},"profils",ybot, ytop);
}

event droplets(t+=1. ; t <= EndTime){

/**
First, we tag the connected bubbles;
  */
	scalar temp[];
	foreach()
		temp[] = f[] >A;
	int nb = tag (temp);
  /**
The `tag()` function identified regions:
Now we see what regions appear at the border of interest.
   */
  	int * bttm = calloc (nb + 1, sizeof(int));;
  	foreach_boundary(bottom)
    		bttm[(int)temp[]] = 1;
	# if _MPI
	  MPI_Allreduce (MPI_IN_PLACE, bttm, nb + 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
	#endif
  /**
The cells with the relevant tag values obtain the original volume
fraction.
   */
	foreach() 
		temp[] = bttm[(int)temp[]] == 1 ? f[] : 0;
	free (bttm);
	boundary ({temp});

	output_ppm (temp, file="bottom-layer.mp4", n = 1024, box = {{-0.5,-0.5},{-0.5 + L0, 0.5}}, map = gray, min = 0., max= 1.); 
  	double area =0.;
  	foreach(reduction(+:area))
	      	if (temp[]==f[]){
			area += temp[]*dv();
	      	}
  	static FILE * fma = fopen("bottom","w");
  	if (t<1.){
    		fprintf (fma, "#[1]t [2]area [3]mixture\n");
  	}
  	fprintf (fma, "%g %g %g\n", t, area, A*L0-area);
  	fflush (fma); 
#if 1
  /** 
  Calculation maximum peak to peak amplitude */
  	static FILE * fa = fopen ("amplitude", "w");
  	if (t<1.){
    		fprintf (fa, "#[1]t [2]Amplitude [3]max [4]min [5]Ra [6]Pr [7]B [8]A\n");
  	}
  	scalar pos[];
	position (temp, pos, {0.,1.});
	stats s = statsf(pos);
	fprintf (fa, "%g %g %.9g %.9g %1.1e %.1f %.2f %.1f\n",
		   t, s.max - s.min, s.max, s.min, Ra, Pr, B, A);
	fflush (fa); 
#endif
}
#endif
/**

## Results

~~~gnuplot
set output "Nu=f(t).png"
set xlabel 'times t'
set ylabel 'Nusselts'
titre=system("awk 'NR==2 {print $2}' data")
set title sprintf("Ra=%s",titre)
plot 'data' u 1:9 w l title 'Nu-top',\
      '' u 1:10 w l title 'Nu-bot'
~~~

![Temperature field.](rbi2dsfdiff/temperature.mp4)
![Animation of the level of refinement.](rbi2dsfdiff/grid.mp4)
![Fraction field.](rbi2dsfdiff/strat.mp4)
![temperature and stratification fields.](rbi2dsfdiff/temp+strat.mp4)
![Bottom layer evolution.](rbi2dsfdiff/bottom-layer.mp4)
*/
