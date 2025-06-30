/**
# Bubble-mediated gas transfer in turbulence
This is the example discussed in section 4 of [Farsoiya et al.,
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
@article{farsoiya2020,
  title = {Bubble mediated gas transfer in turbulence},
  author = {P. K. Farsoiya and S. Popinet and L. Deike},
  journal = {Journal of Fluid Mechanics},
  year = {2020},
}

~~~   
   */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "henry.h"
#include "view.h"

/** We monitor performance statistics and control the maximum runtime. */

#include "navier-stokes/perfs.h"
#include "maxruntime.h"
#include <sys/stat.h>
/** Include density and viscosity ratios.*/
#define RHOR 850.0
#define MUR 25.0
/** Defined domain size is */
#define WIDTH 120.0

 face vector av[];

/** The tracers for gas diffusion */

scalar cSc1[], cSc10[], cSc100[], * stracers = 
  {cSc1, cSc10,  cSc100};

/**
The code takes the level of refinement as optional command-line
argument (as well as an optional maximum runtime). */

int maxlevel = 5;	//only for the compilation
double MAXTIME = 500;	//end time default
int FORCED = 1;		//equals 0 if it is not forced, 1 otherwise
double amp_force = 0.1; //amplitude of the forcing
double SIGMA = 500.0;   // surface tension coeff 
double R0 = 8.0;	//Radius  of the bubble

double SC1 = 1.0;
double SC2 = 1.0;
double SC3 = 1.0;


double visp = 1;
double gravity = 0;

double conc_liq1 = 0;
double conc_gas1 = 1.0;

int main (int argc, char * argv[]) {
  if (argc > 1)
    {
      maxlevel = atoi(argv[1]);
      MAXTIME = atoi(argv[2]);
      SC1 = atof(argv[3]); 
      SC2 = atof(argv[4]); 
      SC3 = atof(argv[5]); 
      
      FORCED = atoi(argv[6]); 		//equals 1 if forced, 0 otherwise
      amp_force = atof(argv[7]);  	//amplitude of the forcing
       //visp for Re = 38, visp=2, Re = 55, visp = 1 and Re = 77, visp = 0.5
      visp = atof(argv[8]);
      SIGMA = atof(argv[9]); 
	gravity = atof(argv[10]);		
    }
  size (WIDTH);
  //~ foreach_dimension()
  periodic (right);
   
  mu1 = visp* 0.01 * sq(WIDTH/(2.0*pi))/(2.0)/2.0; 
  mu2 = mu1/MUR;
  /**
    Use reference density for fluid, and define gas by density ratio.*/
  rho1 = 1.0;
  rho2 = rho1/RHOR;
 
  cSc1.D1 = mu1/rho1/SC1; //diffusivity of gas 1 in liquid
  
  cSc10.D1 = mu1/rho1/SC2; //diffusivity of gas 3 in liquid
 

  cSc100.D1 = mu1/rho1/SC3; //diffusivity of gas 6 in liquid

  cSc1.D2 = cSc1.D1*100.;  //diffusivity of gas 1 in bubble
  
  cSc10.D2 = cSc10.D1*100.;  //diffusivity of gas 3 in bubble
  
  cSc100.D2 = cSc100.D1*100.;  //diffusivity of gas 6 in bubble

  cSc1.alpha=0.3;
  
  cSc10.alpha=0.3;
  
  cSc100.alpha=0.3;

  /**
    And the surface tension.*/

  f.sigma = SIGMA;
  a = av;


#if TREE
  N = 1 << (maxlevel-2);
#else
  N = 1 << maxlevel;
#endif

  /**
    Set tolerance on flow divergence. */
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
   run();
}

/**
## Initial conditions

The initial condition is "ABC" flow. This is a laminar base flow that 
is easy to implement in both Basilisk and a spectral code. */

event init (i = 0) {
  if (!restore (file = "restart") && !restore(file="dump"))
    {
      double waveno = WIDTH/(2.0*pi);
      double amp = WIDTH/(2.0*pi);
      foreach() {

	double xw = x/waveno;
	double yw = y/waveno;
	double zw = z/waveno;
//	u.x[] = 10;
	//u.x[] = amp*( cos(yw) + sin(zw) );
	//u.y[] = amp*( sin(xw) + cos(zw) );
	//u.z[] = amp*( cos(xw) + sin(yw) );
      }
	//~ refine(sq(2.0*R0)-sq(x-WIDTH*0.9) - sq(y-WIDTH/2.0) - sq(z-WIDTH/2.0)>0 && level < maxlevel);
      //~ fraction (f, sq(x-WIDTH*0.9) + sq(y-WIDTH/2.0) + sq(z-WIDTH/2.0) - sq(R0));

    }
  else if(restore (file = "dump")){
    // do nothing
	    }
  else
    {
      refine(sq(2.0*R0)-sq(x-WIDTH/2.0) - sq(y-WIDTH/2.0) - sq(z-WIDTH/2.0)>0 && level < maxlevel);
      fraction (f, sq(x-WIDTH/2.0) + sq(y-WIDTH/2.0) + sq(z-WIDTH/2.0) - sq(R0));
      foreach(){      
	foreach_dimension(){
	  u.x[] = f[]*u.x[];
	}
	
	cSc1[] = conc_liq1*f[] + conc_gas1*(1. - f[]);
	cSc10[] = conc_liq1*f[] + conc_gas1*(1. - f[]);
	cSc100[] = conc_liq1*f[] + conc_gas1*(1. - f[]);

		
      }
	
	  char namedump[80];
  	 sprintf(namedump,"./dat/dump_%g",t);
  	dump(file = namedump);

    }
}

/**
Include accelerations:*/
event acceleration (i++) {

	
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
      
      foreach_face(x){
	  av.x[] += gravity;
	//~ printf("\ngravity = %g\n",av.x[]);
	      }
	
	//~ exit(0);
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
  //~ fprintf (fd, "%g %g %g %g\n",
	   //~ t, vd, ke, 2./3.*ke/mu1*sqrt(15.*mu1/vd));
  fflush (fd);
}

/**
We use adaptivity. */
#if TREE
event adapt (i++) {
  double uemax = 0.2*normf(u.x).avg;
  double tr1emax = 0.2*normf(cSc1).avg;
  double tr3emax = 0.2*normf(cSc10).avg;
  double tr6emax = 0.2*normf(cSc100).avg;

  double femax = 1e-2;

  adapt_wavelet ((scalar *){f, u,cSc1,cSc10,cSc100}, 
		 (double[]){femax, uemax, uemax, uemax,tr1emax,tr3emax,tr6emax},maxlevel);


}
#endif


/**
We output a full snapshot every time unit. */

event snapshot (t=0; t <=MAXTIME; t+=0.01)
{

  char namedump[80];
  sprintf(namedump,"./dat/dump_%g",t);
  //~ dump(file = namedump);
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

	
  double tr1 = 0., tr2 = 0., tr3 = 0.,tr1w=0.,tr2w=0., tr3w = 0.;
  double tr4 = 0., tr5 = 0., tr6 = 0.,tr4w=0.,tr5w=0., tr6w = 0.;
	 
  double vbx = 0.,vby=0., sb=0., vbz = 0.;
  double xb = 0., yb = 0., zb =0. ;
  double area=0;

	
  foreach(
	  reduction(+:tr1) reduction(+:tr2) reduction(+:tr3) 
	  reduction(+:tr1w) reduction(+:tr2w) reduction(+:tr3w)
	  //~ reduction(+:sb) reduction(+:sbi) reduction(+:area) 
	  // nreductmax violated, moved to next iterator
	  ) {
		
    tr1 += cSc1[]*(1. - f[])/(f[]*cSc1.alpha + (1. - f[]))*dv();
    
    tr3 += cSc10[]*(1. - f[])/(f[]*cSc10.alpha + (1. - f[]))*dv();
	
    tr1w+= cSc1[]*f[]*cSc1.alpha/(f[]*cSc1.alpha + (1. - f[]))*dv(); 
    
    tr3w+= cSc10[]*f[]*cSc10.alpha/(f[]*cSc10.alpha + (1. - f[]))*dv(); 
 
  }

  foreach(
	  reduction(+:tr4) reduction(+:tr5) reduction(+:tr6) 
	  reduction(+:tr4w) reduction(+:tr5w) reduction(+:tr6w)
	  ) {
	   
   
    tr6 += cSc100[]*(1. - f[])/(f[]*cSc100.alpha + (1. - f[]))*dv();
	
    
    tr6w += cSc100[]*f[]*cSc100.alpha/(f[]*cSc100.alpha + (1. - f[]))*dv(); 

		 
  }
  foreach(

	  reduction(+:sb) reduction(+:area)
	  reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
	  reduction(+:xb) reduction(+:yb) reduction(+:zb)

	  ) {
			
    double dvb = (1. - f[])*dv();
    sb += dvb;

    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f), p;

      double alpha = plane_alpha (f[], n);
	
      area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    }	

	//location and velocity of the bubble. Not being used now	
    xb += x*dvb;
    yb += y*dvb;
    zb += z*dvb;

    vbx += u.x[]*dvb;
    vby += u.y[]*dvb;
    vbz += u.z[]*dvb;

  }


  static FILE * fp22 = fopen ("gas_transfer_stats.dat", "w");

  if (i == 0)
    fprintf(fp22,"#t tr1 tr2 tr3 tr4 tr5 tr6 tr1w tr2w tr3w tr4w tr5w tr6w sb sbi area vd ke ReL");
	
  //~ fprintf(fp22,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e " 
	  //~ "%e %e %e %e\n",t,tr1,tr2,tr3,tr4,tr5,
	  //~ tr6,tr1w,tr2w,tr3w,tr4w,tr5w,tr6w,sb,statsf(f).sum,area,vd,ke,2./3.*ke/mu1*sqrt(15.*mu1/vd));
	
  fflush (fp22);


}

event movie (t = 0; t +=0.01; t<100)
{
  char namedump[80];
  sprintf(namedump,"Images/snap_%g.ppm",t*100);
 
 view (fov = 44, camera = "iso", ty = .2,
	width = 600, height = 600,  samples = 4);
  clear();
  box();
	squares("p",alpha=60);
  //~ draw_vof("f");
  save (namedump);
}

/**
End the simulation. */

event end (t=MAXTIME) {
  dump ("end");
}

/** ## Results 

~~~pythonplot Sherwood number vs time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

levich = lambda pe: 3.0**0.25*2.*np.sqrt(pe/np.pi);

def sherwood(filename,mu,t0):
    
    Sc1=1; Sc2=2; Sc3=10; Sc4=20; Sc5=50; Sc6=100;
    
    dl1 = 0.913*mu/Sc1;dl2 = 0.913*mu/Sc2;dl3 = 0.913*mu/Sc3; dl4 = 0.913*mu/Sc4;
    
    dl5 = 0.913*mu/Sc5;dl6 = 0.913*mu/Sc6; r0=8; sigma=500;
    
    He1=0.3;He2=0.3;He3=0.3;He4=0.3;He5=0.3;He6=0.3;

    pi=3.1416;
    rho=1;
    t1,tr1,tr2,tr3,tr4,tr5,tr6,tr1w,tr2w,tr3w,tr4w,tr5w,tr6w,vol,sbi,area,vd,ke,ReL = np.loadtxt(filename,delimiter=' ',unpack=True)
    tc=np.power(2.*r0,2/3)*np.power(vd[-1],-1/3);

    urms = np.sqrt(2.*ke/3);
    ti1=[];
    veli=[];
    voli=[];
    areai=[];
    urmsi=[];
    vdi=[];
    kei=[];
    ReLi=[];
    ci1=[];ci2=[];ci3=[];ci4=[];ci5=[];ci6=[];
    co1=[];co2=[];co3=[];co4=[];co5=[];co6=[];
    volw=[];
    trt1=[];trt2=[];trt3=[]; #total moles sanity check
    sh1=[];sh2=[];sh3=[];sh4=[];sh5=[];sh6=[];

    kh=[];
    levich1=[];levich2=[];levich3=[];levich4=[];levich5=[];levich6=[];

    #filename = 'ch_3.dat'
    #th,trh,velh,volh,areah = np.loadtxt(filename,delimiter='\t',unpack=True)
    #for i in range(0, len(t1)): 
       # trt1.extend([tr1[i] + tr1w[i]]); 
       # trt2.extend([tr2[i] + tr2w[i]]); 
       # trt3.extend([tr3[i] + tr3w[i]]); 

    for i in range(0,len(t1)-1):
        ti1.extend([(0.5*(t1[i+1]+t1[i])-t0)/tc]);

        voli.extend([0.5*(vol[i+1]+vol[i])]);
        areai.extend([0.5*(area[i+1]+area[i])]);
        urmsi.extend([0.5*(urms[i+1]+urms[i])]);
        vdi.extend([0.5*(vd[i+1]+vd[i])]);
        kei.extend([0.5*(ke[i+1]+ke[i])]);
        ReLi.extend([0.5*(ReL[i+1]+ReL[i])]);

        ci1.extend([0.5*(tr1[i+1]/vol[i+1]+tr1[i]/vol[i])]);
        ci2.extend([0.5*(tr2[i+1]/vol[i+1]+tr2[i]/vol[i])]);
        ci3.extend([0.5*(tr3[i+1]/vol[i+1]+tr3[i]/vol[i])]);
        ci4.extend([0.5*(tr4[i+1]/vol[i+1]+tr4[i]/vol[i])]);
        ci5.extend([0.5*(tr5[i+1]/vol[i+1]+tr5[i]/vol[i])]);
        ci6.extend([0.5*(tr6[i+1]/vol[i+1]+tr6[i]/vol[i])]);

        co1.extend([0.5*(tr1w[i+1]/sbi[i+1]+tr1w[i]/sbi[i])]);
        co2.extend([0.5*(tr2w[i+1]/sbi[i+1]+tr2w[i]/sbi[i])]);
        co3.extend([0.5*(tr3w[i+1]/sbi[i+1]+tr3w[i]/sbi[i])]);
        co4.extend([0.5*(tr4w[i+1]/sbi[i+1]+tr4w[i]/sbi[i])]);
        co5.extend([0.5*(tr5w[i+1]/sbi[i+1]+tr5w[i]/sbi[i])]);
        co6.extend([0.5*(tr6w[i+1]/sbi[i+1]+tr6w[i]/sbi[i])]);


        sh1.extend([abs((tr1[i+1]-tr1[i])/((t1[i+1]-t1[i]))/areai[i]/(ci1[i]*He1-co1[i]))*2.0*r0/dl1]);
        sh2.extend([abs((tr2[i+1]-tr2[i])/((t1[i+1]-t1[i]))/areai[i]/(ci2[i]*He2-co2[i]))*2.0*r0/dl2]);
        sh3.extend([abs((tr3[i+1]-tr3[i])/((t1[i+1]-t1[i]))/areai[i]/(ci3[i]*He3-co3[i]))*2.0*r0/dl3]);
        sh4.extend([abs((tr4[i+1]-tr4[i])/((t1[i+1]-t1[i]))/areai[i]/(ci4[i]*He4-co4[i]))*2.0*r0/dl4]);
        sh5.extend([abs((tr5[i+1]-tr5[i])/((t1[i+1]-t1[i]))/areai[i]/(ci5[i]*He5-co5[i]))*2.0*r0/dl5]);
        sh6.extend([abs((tr6[i+1]-tr6[i])/((t1[i+1]-t1[i]))/areai[i]/(ci6[i]*He6-co6[i]))*2.0*r0/dl6]);
    
    #calculate sherwood from the farsoiya et. al. equation
    t1 = (t1 - t0)/tc;
    
    
    far_sh1 = levich(urms*2*r0/dl1)
    far_sh2 = levich(urms*2*r0/dl2)
    far_sh3 = levich(urms*2*r0/dl3)
    far_sh4 = levich(urms*2*r0/dl4)
    far_sh5 = levich(urms*2*r0/dl5)
    far_sh6 = levich(urms*2*r0/dl6)
#     print(vd)
    return ti1,urms,sh1,sh2,sh3,sh4,sh5,sh6,t1,far_sh1,far_sh2,far_sh3,far_sh4,far_sh5,far_sh6

ti1,urms,sh1,sh2,sh3,sh4,sh5,sh6,t1,far_sh1,far_sh2,far_sh3,far_sh4,far_sh5,far_sh6 = sherwood('gas_transfer_stats.dat',0.5,116);
plt.semilogy(ti1, sh1,'g',label='Sc=1 old')
# print(urms)
plt.plot(ti1, sh2,'r',label='Sc=2')

plt.plot(ti1, sh3,'m',label='Sc=10')
plt.plot(ti1, sh4,'b',label='Sc=20')

plt.plot(ti1, sh5,'k',label='Sc=50')
plt.plot(ti1, sh6,'y',label='Sc=100')

plt.plot(t1, far_sh1,'g-.')
plt.plot(t1, far_sh2,'r-.')
plt.plot(t1, far_sh3,'m-.')
plt.plot(t1, far_sh4,'b-.')
plt.plot(t1, far_sh5,'k-.')
plt.plot(t1, far_sh6,'y-.')      

plt.ylim(10, 1000)
plt.xlim(0, 1);

plt.xlabel(r'$(t-t_0)/t_c$',fontsize=20)
plt.ylabel(r'$Sh$',rotation=0,fontsize=20)
plt.tight_layout()

plt.savefig('fig7b.png')

~~~
*/
