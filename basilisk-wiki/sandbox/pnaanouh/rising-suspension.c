/**
# Simulation of a rising suspension of drops 

*/ 
/** Necessary header files. no-coalescence.h can be found [here](http://basilisk.fr/sandbox/popinet/no-coalescence.h)*/
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "no-coalescence.h"
#include "view.h"
#include "output.h"
/** cgposition_opt.h  is used to calculate the global parameters for each drop requires [this patch by Antoon van Hooft](http://basilisk.fr/sandbox/Antoonvh/reductions?raw) for array reduction with MPI */
#include "cgposition_opt.h"

/** We define the characteristic parameters of the flow where Bo is the Bond/Eötvös number, Ga is the Galileo number equivalent to $\sqrt{Ar}$ the Archimedes number, rhor and rmu are respectively the density and viscosity ratios between the continuous and dispersed phases, and PHI is the volume fraction  */
#define Bo 0.5
#define Ga 55.
#define rhor 1.1668599
#define rmu 0.4237805
#define PHI 0.15
/** We define N to be the number of unit cells per dimension of the domain*/
#define N 5
/** D is the diameter of the drops is set to one for convenience. The Domain size Ls and the unit cell size DP are calculated as functions of N, D, and PHI. RND is the maximum displacement of the drop from the center of the unit cell at initialization */
#define D 1.
#define Ls N*D*sqrt(pi/(4.*PHI))
#define DP Ls/N
#define RND 0.49*(DP-D)

/** We set the gravitational acceleration and the density of the dispersed phase to one for convenience */
#define grav 1.
#define r1 1.
/** We define the maximum simulation time and mesh refinement level*/
#define TMAX 2000
#define LEVEL  9

/** Function used to initialise the drops on a square grid so the maximum possible volume fraction $\phi=\frac{\pi}{4}\approx0.785$*/
double geometry( double x, double y, double z);

double geometry( double x, double y, double z){
  double inout=-1.,circ=0.;
  srand(3458829);
  for(double i=(-Ls/2.+DP/2.); i<= Ls/2; i+=DP){
    for( double j=(-Ls/2.+DP/2.); j<= Ls/2.; j+=DP){     
      int A=rand();
      int B=rand();
      circ= -sq(x-i+RND*(A*2./RAND_MAX-1))-sq(y-j+RND*(B*2./RAND_MAX-1))+sq(D/2);
      if (circ>=0){
	inout=circ;
      }									    
    }
  }
  return(inout);
}

/** Event to recover the global parameters of each drop works by tagging each drop and then outputs to stderr the iteration, time, position of the CG, velocity, and tag value for the drop*/
scalar tag_level[],tg_l[];
int nbub;


event track_bub(t+=0.1){
  foreach()
    tg_l[]=0;
  for(scalar s in interfaces){
    foreach()
      tag_level[]=s[]>1e-4;
    nbub=tag(tag_level);
    
    cg poscg[nbub];
    cg_bub(s,tag_level,nbub,poscg);
    for(int j=0;j<nbub;j++){
      if(poscg[j].vol>1e-1){
	fprintf(stderr,"%d %g %g %g  %g %g %g %d\n",i,t,poscg[j].x,poscg[j].y,poscg[j].vx,poscg[j].vy,poscg[j].vol,j+1);
      }
      else
        foreach(){ //Fragment suppression 
          if(tag_level[]==j+1)
            s[]=0;
        }
    }
    foreach(){
      tag_level[]*=s[];
      tg_l[]+=tag_level[];
    }
  }
} 

/** The average(gl) of the flow properties over the entire computational domain is computed and outputed to stdout. */
/*
The outputs are: 
* i: iterations
* t: time
* vx/yd: velocity of the dispersed phase
* vx/yf: velocity of the continous phase
* rhox/y: total momentum 
* Vx/y: velocity of the computational domain

*/

event average_stats(t+=0.1){
  double vxd=0,vyd=0,vxf=0,vyf=0,rhovx=0,rhovy=0,mt=0,vtd=0,vtf=0;
  if (i==0)
    fprintf(stdout,"i t N_Vof vxdrops vydrops vxfluid vyfluid rhovx rhovy Vx Vy\n");
  foreach(reduction(+:vxd) reduction(+:vyd) reduction(+:vxf) reduction(+:vyf) reduction(+:rhovx) reduction(+:rhovy) reduction(+:mt) reduction(+:vtd) reduction(+:vtf) ){
    double massl=dv()*(f[]*(rho1-rho2)+rho2);
    vxd+=dv()*f[]*u.x[];
    vyd+=dv()*f[]*u.y[];
    vxf+=dv()*(1-f[])*u.x[];
    vyf+=dv()*(1-f[])*u.y[];
    rhovx+=massl*u.x[];
    rhovy+=massl*u.y[];
    vtd+=dv()*f[];
    vtf+=dv()*(1-f[]);
    mt+=massl;
  }
  double Vx=0,Vy=0;
  foreach(reduction(+:Vx) reduction(+:Vy)){
    Vx+=dv()*u.x[];
    Vy+=dv()*u.y[];
  }
  vxd=vxd/vtd;
  vyd=vyd/vtd;
  vxf=vxf/vtf;
  vyf=vyf/vtf;
  rhovx=rhovx/mt;
  rhovy=rhovy/mt;
  Vx=Vx/(Ls*Ls);
  Vy=Vy/(Ls*Ls);
  int nvof=0;
  for(scalar s in interfaces ){
    reduction(+:nvof);
    nvof+=1;
  }
  fprintf(stdout,"%d %g %d %g %g %g %g %g %g %g %g \n",i,t,nvof,vxd,vyd,vxf,vyf,rhovx,rhovy,Vx,Vy);
}


/** The flow parameters are calculated based on the previously defined non-dimensional parameters. The tolerance is reduced to 1e-4 and the boundaries are set to be periodic*/
int main(){
  size(Ls); 
       init_grid(1 << (LEVEL));
  origin(-L0/2.,-L0/2.,-L0/2.);
  rho1=r1;
  mu1=sqrt(fabs(1-rhor)*rhor*grav*(D*D*D))*rho1/(Ga*rmu);
  rho2=rho1*rhor;
  mu2=mu1*rmu;
  f.sigma=fabs(1-rhor)*rho1*grav*D*D/Bo;

  TOLERANCE=1e-4;
  
  foreach_dimension()
    periodic(right);
  run();
}
/**Initializes the domain with zero velocity and outputs the flow parameters*/
event init (t=0){
  fraction(f, geometry(x,y,z));
  foreach()
    u.x[]=u.y[]=0.;  
  fprintf(stderr,"rho1 %g rho2 %g mu1 %g m2 %g sigma %g Ls %g L0 %g \n",rho1,rho2,mu1,mu2,f.sigma,Ls,L0);
}


/** Overloads the acceleration event to apply the effect of gravity. Due to the periodic boundary conditions the acceleration needs to be reduced by $\frac{\rho_{av}}{\rho}g$ rho is not available directly as a face centered vector so the stagger of f[] is adjusted to calculate it */
event acceleration (i++) {
  double rhoav=avg_rho();
  face vector av = a;
  foreach_face(y)
    av.y[] -= (1-rhoav/((f[]+f[0,-1])/2*(rho1-rho2)+rho2))*grav; 
}
/** Outputs videos of the velocity, Pressure, vorticity, and tag fields*/
event movies(t=0; t<=TMAX; t+=0.04)
{
  
  box();
  draw_vof("f", lw = 2.);
  squares("u.y");
  save("uy.mp4");

  clear();

  box();
  draw_vof("f", lw = 2.);
  squares("u.x");
  save("ux.mp4");

  clear();

  box();
  draw_vof("f", lw = 2.);
  squares("p");
  save("P.mp4");

  clear();
  
  scalar omega[];
  vorticity (u, omega);

  box();
  draw_vof("f", lw = 2.);
  squares("omega");
  save("vort.mp4");

  clear();
  
  box();
  draw_vof("f", lw = 2.);
  squares(color="tg_l",map=randomap,min=0,max=25);
  save("tag.mp4");
  clear();
}

event stop(t=TMAX){
  return 1;
}
 
