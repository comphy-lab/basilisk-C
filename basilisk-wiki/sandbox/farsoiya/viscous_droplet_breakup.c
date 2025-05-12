/**
# Droplet breakage in turbulent flow 
Role of viscous droplet breakup for the [Farsoiya et. al. 2023](#farsoiya2023role). 

*/
//~ #define MIN_LEVEL 5
//~ #define LEVEL 10
int MAX_LEVEL = 8;	//only for the compilation
//~ #define dR_refine (2.*L0/(1 << LEVEL))

#define F_ERR 1e-10
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

/** We monitor performance statistics and control the maximum runtime. */

#include "navier-stokes/perfs.h"
#include "maxruntime.h"
#include <sys/stat.h>
/** Include density and viscosity ratios.*/
#define RHOR 1.
//#define MUR 0.004
/** Defined domain size is */
#define WIDTH 120.0

 


/**
The code takes the level of refinement as optional command-line
argument (as well as an optional maximum runtime). */

double MAXTIME = 500;	//end time default
int FORCED = 1;		//equals 0 if it is not forced, 1 otherwise
double amp_force = 0.1; //amplitude of the forcing
double SIGMA = 500.0;   // surface tension coeff 
double R0 = 8.0;	//Radius  of the bubble

double SC1 = 1.0;
double mur = 1.;
double We =1.;


double visp = 1.;
double gravity = 0;
double tst = 1;
double conc_liq1 = 0;
//double conc_gas1 = 1.0;

int main (int argc, char * argv[]) {
  if (argc > 1)
    {
      MAX_LEVEL = atoi(argv[1]);
      MAXTIME = atoi(argv[2]);
      SC1 = atof(argv[3]); 
 

      FORCED = atoi(argv[4]); 		//equals 1 if forced, 0 otherwise
      amp_force = atof(argv[5]);  	//amplitude of the forcing
					  //visp for Re = 38, visp=2, Re = 55, visp = 1 and Re = 77, visp = 0.5
							  visp = atof(argv[6]);
      We = atof(argv[7]); 
      mur = atof(argv[8]);
      tst = atof(argv[9]);		
    }
  size (WIDTH);
  foreach_dimension()
    periodic (right);
  /**
    Use reference density for fluid, and define gas by density ratio.*/
  double epsilon = 10.;
  rho1 = 1.0;
  rho2 = rho1/RHOR;  
  SIGMA = 2.*rho1*pow(epsilon,2./3)*pow(2.*R0,5./3)/We;
  mu1 = visp* 0.01 * sq(WIDTH/(2.0*pi))/(2.0)/2.0;  
  mu2 = mu1*mur; 

  /**
    And the surface tension.*/

  f.sigma = SIGMA;
  


#if TREE
  N = 1 << (MAX_LEVEL-2);
#else
  N = 1 << MAX_LEVEL;
#endif

  /**
    Set tolerance on flow divergence. */
    TOLERANCE = 1e-4;
    NITERMAX = 20;
    
    struct stat st = {0};
    char name[80];
    sprintf (name, "dat");
    char newname[80];
    sprintf (name, "dat");

    if (stat(name, &st) == -1)
      {
	mkdir(name, 0700);
      }
    sprintf (name, "Images");

    if (stat(name, &st) == -1)
      {
	mkdir(name, 0700);
      }


    /**
      Preserve old stats restoring the simulation with a timestamp suffixed. */		
      sprintf (name, "stats.dat");
	
    if (stat(name, &st) == 0)
      {
	sprintf (newname, "stats-%ld.dat",time(0));
	rename(name, newname);
      }
    sprintf (name, "perfs");
	
    if (stat(name, &st) == 0)
      {
	sprintf (newname, "perfs-%ld",time(0));
	rename(name, newname);
      }

    sprintf (name, "gas_transfer_stats.dat");
	
    if (stat(name, &st) == 0)
      {
	sprintf (newname, "gas_transfer_stats-%ld.dat",time(0));
	rename(name, newname);
      }	
    sprintf (name, "bubbles.dat");
	
    if (stat(name, &st) == 0)
      {
	sprintf (newname, "bubbles-%ld.dat",time(0));
	rename(name, newname);
      }
    sprintf (name, "droplets.dat");
	
    if (stat(name, &st) == 0)
      {
	sprintf (newname, "droplets-%ld.dat",time(0));
	rename(name, newname);
      }
 
    run();
}

/**
## Initial conditions

The initial condition is "ABC" flow. This is a laminar base flow that 
is easy to implement in both Basilisk and a spectral code. */
#include <dirent.h>
event init (i = 0) {

  if (!restore (file = "restart") && !restore(file="dump"))
    {
      double waveno = WIDTH/(2.0*pi);
      double amp = WIDTH/(2.0*pi);
      foreach() {
	f[] = 1.;
	double xw = x/waveno;
	double yw = y/waveno;
	double zw = z/waveno;

	u.x[] = amp*( cos(yw) + sin(zw) );
	u.y[] = amp*( sin(xw) + cos(zw) );
	u.z[] = amp*( cos(xw) + sin(yw) );
      }

    }
  else if(restore (file = "dump")){
    // do nothing
	    }
  else
    {
      refine(sq(2.0*R0)-sq(x-WIDTH/2.0) - sq(y-WIDTH/2.0) - sq(z-WIDTH/2.0)>0 && level < MAX_LEVEL );
      fraction (f, sq(x-WIDTH/2.0) + sq(y-WIDTH/2.0) + sq(z-WIDTH/2.0) - sq(R0));
      foreach(){      
	foreach_dimension(){
	  u.x[] = f[]*u.x[];
	}
	
		
      }
	


    }
}

/**
Include accelerations:*/
event acceleration (i++) {

  face vector av = a;

  /** ## Linear forcing term
    We compute the average velocity and add the corresponding linear
    forcing term. */
    if (FORCED)
      {
	coord ubar;
	foreach_dimension() {
	  stats s = statsf(u.x);
	  ubar.x = s.sum/s.volume;
	}

	foreach_face()
	  av.x[] += f[]*amp_force*((u.x[] + u.x[-1])/2. - ubar.x);

      }

  boundary((scalar *){av});
}


/**
## Outputs

We log the evolution of viscous dissipation, kinetic energy, and
  microscale Reynolds number. */

  event logfile (i++) {
    //~ printf("\n%d\t%g\t%g",i,t,dt);
    //~ fflush(stdout);
    coord ubar;
    foreach_dimension() {
      stats s = statsf(u.x);
      ubar.x = s.sum/s.volume;
    }

    double ke = 0., vd = 0., vol = 0.;
    foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
      vol += dv();
      foreach_dimension() {
	// mean fluctuating kinetic energy
	ke += dv()*sq(u.x[] - ubar.x);
	// viscous dissipation
	vd += dv()*(sq(u.x[1] - u.x[-1]) +
		   sq(u.x[0,1] - u.x[0,-1]) +
		   sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
      }
    }
    ke /= 2.*vol;
    vd *= mu1/vol;
    static FILE * fd = fopen("stats.dat","w");
    if (i == 0) {
      //~ fprintf (ferr, "t dissipation energy Reynolds\n"); 
      fprintf (fd, "t dissipation energy Reynolds\n");
    }
    //~ fprintf (ferr, "%g %g %g %g\n",
		 //~ t, vd, ke, 2./3.*ke/mu1*sqrt(15.*mu1/vd));
    fprintf (fd, "%g %g %g %g\n",
	     t, vd, ke, 2./3.*ke/mu1*sqrt(15.*mu1/vd));
    fflush (fd);
  }

/**
We use adaptivity. */
#if TREE
event adapt (i++) {
  double uemax = 0.2*normf(u.x).avg;
  // double tr1emax = 0.2*normf(c1).avg;


  double femax = 1e-2;

  adapt_wavelet ((scalar *){f, u}, 
		(double[]){femax, uemax, uemax, uemax},MAX_LEVEL);


}
#endif


/**
We output a full snapshot every time unit. */
double endtimemark = 300;

event snapshot (t=0; t <=MAXTIME; t+=1)
{

  char namedump[80];
  sprintf(namedump,"./dat/dump_%g",t);
  dump(file = "dump");
}

/**
We also want to count the drops and bubbles in the flow. */
int endflag = 0;
event countDropsBubble(t+=0.01)
{
  scalar m1[]; //droplets
  scalar m2[]; //bubbles
  foreach(){
    m1[] = f[] > 0.05; //i.e. set m true if f[] is close to unity (droplets)
    m2[] = f[] < 0.99; //m true if f[] close to zero (bubbles)
				      }
  int n1 = tag(m1);
  int n2 = tag(m2);
  /**
    Having counted the bubbles, now we find their size. This example
    is similar to the jet atomization problem. We are interested in
    the volumes and positions of each droplet/bubble.*/
    double v1[n1]; //droplet
  coord b1[n1];  //droplet
  double v2[n2]; //bubble
  coord b2[n2];  //bubble
  double ar1[n1];
  double ar2[n2];
  /**
    We initialize: */
    for (int j=0; j<n1; j++)
      {
	v1[j] = ar1[j] = b1[j].x = b1[j].y = b1[j].z = 0.0;
      }
  for (int j=0; j<n2; j++)
    {
      v2[j] = ar2[j] = b2[j].x = b2[j].y = b2[j].z = 0.0;
    }
  /**
    We proceed with calculation. */
    foreach_leaf() //droplets
    {
      if (m1[] > 0) {
	int j = m1[] - 1;
	v1[j] += dv()*f[]; //increment the volume of the droplet
	if (f[] > 1e-6 && f[] < 1. - 1e-6) {
	  coord n = interface_normal (point, f), p;

	  double alpha = plane_alpha (f[], n);
	
	  ar1[j] += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
	}	

 
     
	coord p = {x,y,z};
	foreach_dimension()
	  b1[j].x += dv()*f[]*p.x;
      }
    }
  foreach_leaf() //bubbles
    {
      if (m2[] > 0) {
	int j = m2[] - 1;
	v2[j] += dv()*(1.0-f[]);
	if (f[] > 1e-6 && f[] < 1. - 1e-6) {
	  coord n = interface_normal (point, f), p;

	  double alpha = plane_alpha (f[], n);
	
	  ar2[j] += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
	}	

	coord p = {x,y,z};
	foreach_dimension()
	  b2[j].x += dv()*(1.0-f[])*p.x;
      }
    }
  /**
    Reduce for MPI. */
#if _MPI
		      MPI_Allreduce (MPI_IN_PLACE, v1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, ar1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (MPI_IN_PLACE, b1, 3*n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, v2, n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, ar2, n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (MPI_IN_PLACE, b2, 3*n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  /**
    Output the volume and position of each droplet to file. */
    static FILE * fdrop = fopen("droplets.dat","w");
  static FILE * fbubb = fopen("bubbles.dat","w");
  for (int j=0; j<n1; j++)
    {
      fprintf (fdrop, "%d %g %d %g %g %g %g %g\n", i, t,
	       j, v1[j], b1[j].x/v1[j], b1[j].y/v1[j],b1[j].z/v1[j],ar1[j]);
    }
  for (int j=0; j<n2; j++)
    {
      fprintf (fbubb, "%d %g %d %g %g %g %g %g\n", i, t,
	       j, v2[j], b2[j].x/v2[j], b2[j].y/v2[j], b2[j].z/v2[j],ar2[j]);
    }
  double temp;
  // Simulation end criteria
    if (n2 > 1){   //if more than one droplet exists
	//sort droplet volume ascending
	for (int k=0; k<n2; k++){
	  for (int j=0; j<(n2-k); j++){
	    if(v2[j] > v2[j+1]){
	      temp = *(v2 + j);
	      *(v2 + j) = *(v2 + j + 1);
	      *(v2 + j + 1) = temp;
	    }
	  }
	}
      double delta = L0/(pow(2,MAX_LEVEL));
      // Thanks to Alienor who pointed out that looking the size of smallest droplet can result in judging ficticious breakup
	// //hence we will now look for second largest droplet which should be more than 4 delta cube
  if (endflag == 0 && v2[n2-2] > 4.*pow(delta,3)){ //if first breakup and smallest droplet has volume more than 4 times cell volume
  endflag  = 1; //declare fist breakup
  endtimemark = t; //not the time for first breakup
  dump(file = "breakup"); //save the breakup dump file
  }
 else if (endflag == 1){
  if (t > (endtimemark + 30)){ // end simulation after 10tc
  printf("End time after breakup"); 
  dump(file = "end");
  exit(0);
  }
 }
     
 }   	

}




scalar omg[];
void vorticity3D (vector u, scalar omg) {
  foreach() {
    omg[] = 0.;
    foreach_dimension()
      omg[] += sq((u.y[1] - u.y[-1] -
		   u.x[0,1] + u.x[0,-1])/(2*Delta));

    omg[] = sqrt(omg[]);
  }
  boundary ((scalar*){omg});
}

event movie (t+=0.1)
{
  char namedump[200];
  sprintf(namedump,"Images/temp.ppm");
  vorticity3D(u, omg);

  view (fov = 44, camera = "iso", ty = .2,
	width = 1200, height = 1200,  samples = 4);
  clear();
  //box();
  draw_vof("f");
  squares ("omg", spread = 5, linear = true, map = jet, alpha = 1, n = {1, 0 , 0});
  squares ("omg", spread = 5, linear = true, map = jet, alpha = 1, n = {0, 1 , 0});
  squares ("omg", spread = 5, linear = true, map = jet, alpha = 1, n = {0, 0 , 1});
 
  
  save (namedump);
  sprintf(namedump,"convert Images/temp.ppm Images/snap-%g.png",t*100);
  system(namedump);

}

/**
End the simulation. */

event end (t=MAXTIME) {
  dump ("end");
}
/**
## References
~~~bib
@article{farsoiya2023role,
  title={Role of viscosity in turbulent drop break-up},
  author={Farsoiya, Palas Kumar and Liu, Zehua and Daiss, Andreas and Fox, Rodney O and Deike, Luc},
  journal={Journal of Fluid Mechanics},
  volume={972},
  pages={A11},
  year={2023},
  publisher={Cambridge University Press}
}
~~~   
*/