/**
# Bubble rising in turbulence
This is the code used for studying the rising speed of bubble in turbulence [Liu et al.,
2024](#liu2024). . The code is based on the example discussed in section 4 of [Farsoiya et al.,
2020](#farsoiya2020).
Given an initial fully developed turbulent flow, as produced by the example case
  isotropicAdaptive.c based on the case investigated by [Rosales & Meneveau, 2005](/src/references.bib#rosales2005)).

 This simulation is in two parts. First, an equivalent version of isotropic.c (dubbed isotropicAdaptive.c) is run, converted to have consistent length scales with the bubble case but without a bubble present in the flow. Second, the simulation is dumped and restarted, and the bubble is initialized in the flow. 

   Both parts are incorporated into this one problem file using preprocessor commands.

   Include a way to restart the files. The filename "restart" is the dump file from the precursor simulation. 
   To continue the simulation by restoration should have a file dump file called "dump".

   Write several dump files to avoid problems when you want to restart your simulations. Initially flow inside the bubble put to 0. No forcing inside the bubble. 

## References

~~~bib
@article{liu2024,
  title = {Direct numerical simulation of bubble rising
in turbulence},
  author = {Z. Liu and P. K. Farsoiya and S. Perrard and L. Deike},
  journal = {Journal of Fluid Mechanics},
  year = {2024},
  note = {submitted}
}
@article{farsoiya2020,
  title = {Bubble mediated gas transfer in turbulence},
  author = {P. K. Farsoiya and S. Popinet and L. Deike},
  journal = {Journal of Fluid Mechanics},
  year = {2020},
  note = {submitted}
}

~~~   
   */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
#include "lambda2.h"

/** We monitor performance statistics and control the maximum runtime. */
#include "navier-stokes/perfs.h"
#include "maxruntime.h"
#include <sys/stat.h>

/** Include density and viscosity ratios.*/
#define RHOR 850.0
#define MUR 25.0

/** Defined domain size is */
#define WIDTH 120.0

/**
The code takes the level of refinement as optional command-line
argument (as well as an optional maximum runtime). */

int maxlevel = 9;	//only for the compilation
double MAXTIME = 500;	//end time default
int FORCED = 1;		//equals 0 if it is not forced, 1 otherwise
double amp_force = 0.1; //amplitude of the forcing
double SIGMA = 500.0;   // surface tension coeff 
double R0 = 8.0;	//Radius  of the bubble

double visp = 1.;
double gravity = 0;

int main (int argc, char * argv[]) {
  if (argc > 1)
    {
	maxlevel = atoi(argv[1]);
	MAXTIME = atoi(argv[2]);
	FORCED = atoi(argv[3]); 		//equals 1 if forced, 0 otherwise
	amp_force = atof(argv[4]);  	//amplitude of the forcing, we keep it 0.1
	visp = atof(argv[5]); 		//visp for Re = 38, visp=2, Re = 55, visp = 1 and Re = 77, visp = 0.5
	SIGMA = pow(5.78,2)*2.*R0/atof(argv[6]);  // sigma = f(We) find the def of We from paper
	gravity = pow(5.78/atof(argv[7]),2)/R0/2.; // g = f(Fr) find the def from paper, dan Ruth			
    }

  size (WIDTH);
  foreach_dimension()
  periodic (right);
   
  mu1 = visp* 0.01 * sq(WIDTH/(2.0*pi))/(2.0)/2.0; 
  mu2 = mu1/MUR;
  /** Use reference density for fluid, and define gas by density ratio.*/
  rho1 = 1.0;
  rho2 = rho1/RHOR;

  /** And the surface tension.*/
  f.sigma = SIGMA;
  
  #if TREE
    N = 1 << (maxlevel - 2);
  #else
    N = 1 << maxlevel;
  #endif

  /** Set tolerance on flow divergence. */
  TOLERANCE = 1e-4;

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


  /** Preserve old stats restoring the simulation with a timestamp suffixed. */		
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

  sprintf (name, "RB_stats.dat");
	
  if (stat(name, &st) == 0)
    {
      sprintf (newname, "RB_stats-%ld.dat",time(0));
      rename(name, newname);
    }		
   run();
}

/**
## Initial conditions

The initial condition is "ABC" flow. This is a laminar base flow that 
is easy to implement in both Basilisk and a spectral code. */

event init (i = 0) {
  char name[80];
  sprintf(name,"restart_mu_%g",visp);

  if (!restore (file = name) && !restore(file="dump") && FORCED == 1)
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
    // do nothing just for extending the simulation regardless of precursor or with bubble
	    }
  else
    {
      refine(sq(2.0*R0)-sq(x-WIDTH*0.1) - sq(y-WIDTH/2.0) - sq(z-WIDTH/2.0)>0 && level < maxlevel);
      fraction (f, sq(x-WIDTH*0.1) + sq(y-WIDTH/2.0) + sq(z-WIDTH/2.0) - sq(R0));
      foreach(){      
	foreach_dimension(){
	  u.x[] = f[]*u.x[];
	}		
      }
	
  	char namedump[80];
 	sprintf(namedump,"./dat/dump_%g",t);
  	dump(file = namedump);
    }
}

/**
Include accelerations:*/
event acceleration (i++) {
 face vector av = a;

	foreach_face(x)
	  av.x[] = gravity*(-1. + alphav.x[]*rho1);


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
    fprintf (fd, "t dissipation energy Reynolds\n");
  }
  fprintf (fd, "%g %g %g %g\n",
	   t, vd, ke, 2./3.*ke/mu1*sqrt(15.*mu1/vd));
  fflush (fd);
}

/**
We use adaptivity. */
#if TREE
event adapt (i++) {
  double uemax = 0.2*normf(u.x).avg;
  double femax = 1e-2;

  adapt_wavelet ((scalar *){f, u}, 
		 (double[]){femax, uemax, uemax, uemax},maxlevel);


}
#endif


/**
We output a full snapshot every time unit. */

event snapshot (t=0; t <=MAXTIME; t+=1)
{

  char namedump[80];
  sprintf(namedump,"./dat/dump_%g",t);
  dump(file = namedump);
}

event extract (t+=0.01){

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
 
  double vbx = 0.,vby=0., sb=0., vbz = 0.;
  double xb = 0., yb = 0., zb =0. ;
  double xb1 = 0., yb1 = 0., zb1 =0. ;
  double xb2 = 0., yb2 = 0., zb2 =0. ;
  double sbcx = 0., sbcy = 0., sbcz = 0.;
  double area=0;


  foreach(
	  reduction(+:sb) reduction(+:area)
	  reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
	  reduction(+:xb) reduction(+:yb) reduction(+:zb)
  	  reduction(+:xb1) reduction(+:yb1) reduction(+:zb1)
  	  reduction(+:xb2) reduction(+:yb2) reduction(+:zb2)
  	  reduction(+:sbcx) reduction(+:sbcy) reduction(+:sbcz)
	  ) {
			
    double dvb = (1. - f[])*dv();
    sb += dvb;
    sbcx += dvb * (x > 0.25 * WIDTH && x < 0.75 * WIDTH);
    sbcy += dvb * (y > 0.25 * WIDTH && y < 0.75 * WIDTH);
    sbcz += dvb * (z > 0.25 * WIDTH && z < 0.75 * WIDTH);
    

    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f), p;

      double alpha = plane_alpha (f[], n);
	
      area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    }	

    //location and velocity of the bubble. Not being used now	
    xb1 += x*dvb;
    yb1 += y*dvb;
    zb1 += z*dvb;
    
    xb2 += (x + (x < 0.5 * WIDTH) * WIDTH) * dvb;
    yb2 += (y + (y < 0.5 * WIDTH) * WIDTH) * dvb;
    zb2 += (z + (z < 0.5 * WIDTH) * WIDTH) * dvb;
    
    vbx += u.x[]*dvb;
    vby += u.y[]*dvb;
    vbz += u.z[]*dvb;

  }
  
  xb = (sbcx >= 0.5 * sb) * (xb1/sb) + (sbcx < 0.5 * sb) * (fmod(xb2/sb, WIDTH));
  yb = (sbcy >= 0.5 * sb) * (yb1/sb) + (sbcy < 0.5 * sb) * (fmod(yb2/sb, WIDTH));
  zb = (sbcz >= 0.5 * sb) * (zb1/sb) + (sbcz < 0.5 * sb) * (fmod(zb2/sb, WIDTH));
  

  double vsx1 = 0., vsy1 = 0., vsz1 =0., sb1 =0. ;
  double vsx2 = 0., vsy2 = 0., vsz2 =0., sb2 =0. ;
  double vsx3 = 0., vsy3 = 0., vsz3 =0., sb3 =0. ;
  double vsx4 = 0., vsy4 = 0., vsz4 =0., sb4 =0. ;
  double vsx5 = 0., vsy5 = 0., vsz5 =0., sb5 =0. ;
  double vsx6 = 0., vsy6 = 0., vsz6 =0., sb6 =0. ;
  
  foreach(
  	reduction(+:sb1) reduction(+:sb2) reduction(+:sb3)
  	reduction(+:sb4) reduction(+:sb5) reduction(+:sb6)
	reduction(+:vsx1) reduction(+:vsy1) reduction(+:vsz1)
	reduction(+:vsx2) reduction(+:vsy2) reduction(+:vsz2)
	reduction(+:vsx3) reduction(+:vsy3) reduction(+:vsz3)
        reduction(+:vsx4) reduction(+:vsy4) reduction(+:vsz4)
	reduction(+:vsx5) reduction(+:vsy5) reduction(+:vsz5)
	reduction(+:vsx6) reduction(+:vsy6) reduction(+:vsz6)
  ){

    double xb0 = xb, yb0 = yb, zb0= zb;	
    double disx = min(abs(x-xb0),min(abs(x+WIDTH-xb0),abs(x-WIDTH-xb0)));
    double disy = min(abs(y-yb0),min(abs(y+WIDTH-yb0),abs(y-WIDTH-yb0)));
    double disz = min(abs(z-zb0),min(abs(z+WIDTH-zb0),abs(z-WIDTH-zb0)));
    
    double dis = (sq(disx) + sq(disy) + sq(disz)) - sq(R0);
    
    double SR1 = 1.1 * R0;
    double SR2 = 1.2 * R0;
    double SR3 = 1.5 * R0;
    double SR4 = 2.0 * R0;
    double SR5 = 3.0 * R0;
    double SR6 = 5.0 * R0;
    
    double dis1 = sq(SR1) - sq(R0);
    double dis2 = sq(SR2) - sq(R0);
    double dis3 = sq(SR3) - sq(R0);
    double dis4 = sq(SR4) - sq(R0);
    double dis5 = sq(SR5) - sq(R0);
    double dis6 = sq(SR6) - sq(R0);
    
    double dvf = f[] * dv();
    
    if(dis > 0 && dis < dis1){
    	vsx1 += u.x[] * dvf;
    	vsy1 += u.y[] * dvf;
    	vsz1 += u.z[] * dvf;
    	sb1 += dvf;
    }
    
    if(dis > 0 && dis < dis2){
    	vsx2 += u.x[] * dvf;
    	vsy2 += u.y[] * dvf;
    	vsz2 += u.z[] * dvf;
    	sb2 += dvf;
    }
    
    if(dis > 0 && dis < dis3){
    	vsx3 += u.x[] * dvf;
    	vsy3 += u.y[] * dvf;
    	vsz3 += u.z[] * dvf;
    	sb3 += dvf;
    }
    
    if(dis > 0 && dis < dis4){
    	vsx4 += u.x[] * dvf;
    	vsy4 += u.y[] * dvf;
    	vsz4 += u.z[] * dvf;
    	sb4 += dvf;
    }
    
    if(dis > 0 && dis < dis5){
    	vsx5 += u.x[] * dvf;
    	vsy5 += u.y[] * dvf;
    	vsz5 += u.z[] * dvf;
    	sb5 += dvf;
    }
    
    if(dis > 0 && dis < dis6){
    	vsx6 += u.x[] * dvf;
    	vsy6 += u.y[] * dvf;
    	vsz6 += u.z[] * dvf;
    	sb6 += dvf;
    }
  
  }


  static FILE * fp22 = fopen ("RB_stats.dat", "w");
	
  fprintf(fp22,"%.12e %.12e %.12e " 
	  "%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",t,sb,statsf(f).sum,area,vd,ke,2./3.*ke/mu1*sqrt(15.*mu1/vd),xb,yb,zb,vbx/sb,vby/sb,vbz/sb,vsx1/sb1,vsy1/sb1,vsz1/sb1,vsx2/sb2,vsy2/sb2,vsz2/sb2,vsx3/sb3,vsy3/sb3,vsz3/sb3,vsx4/sb4,vsy4/sb4,vsz4/sb4,vsx5/sb5,vsy5/sb5,vsz5/sb5,vsx6/sb6,vsy6/sb6,vsz6/sb6);
	
  fflush (fp22);


}


event movie (t +=0.1)
{
  char namedump[80];
  sprintf(namedump,"Images/snap_%g.ppm",t*100);
 
  view (fov = 44, camera = "iso", ty = .2,
	width = 600, height = 600,  samples = 4);
  clear();
  box();
  draw_vof("f");
  save (namedump);
}

/**
End the simulation. */

event end (t=MAXTIME) {
  dump ("end");
}
