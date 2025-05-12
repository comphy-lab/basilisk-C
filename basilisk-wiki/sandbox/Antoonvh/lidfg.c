/**
# The Lid-Driven Cavitty in 2D
We study the flow in a lid-driven square cavity at different reynolds numbers. The set-up consists of a square ($L \times L$) no-slip box with a top lid that forces a flow by moving in the positive x-direction with velocity $U$. Given a fluids viscosity ($\nu$), we can define a Reynoldsnumber: 

$$Re=\frac{UL}{\nu}.$$

Expessed in advection time scales we run the simulation until $t_{end}=20 \times L/U$.

This set up is the fixed-grid equivalent of the [adaptive grid example set up](lid.c).
*/
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

scalar omg[],lev[],U[];
int maxlevel,m;
double vis,e,ens,vit,vir,vil,vib;

int main(){
  L0=1.;
  X0=Y0=0.;
  /**
  We run the simulation 3 times with an increasing Reynoldsnumber and maximum grid resolution.
  */
  for (m=0;m < 3;m++){
    run();
  }
}
/**
The Boundary conditions are chosen to be a no-slip box with a top lid that drived the flow. 
*/
u.t[top]=dirichlet(1.);
u.t[bottom]=dirichlet(0.);
u.n[left]=dirichlet(0.);
u.n[right]=dirichlet(0.);
u.t[left]=dirichlet(0.);
u.t[right]=dirichlet(0.);
u.n[top]=dirichlet(0.);
u.n[bottom]=dirichlet(0.);

event init(t=0){
  /**
  The Reynoldsnumber is directly controlled by the viscosity of the fluid $\nu$.The Reynolds number and maximum grid resolution are doubled each run iteration, starting with $Re=125$ and a  effective resolution of $64\times64$ grid points.  
  */
  vis=1./(125.*pow(2.,(double)m));
  maxlevel=6+m;
  /**
  The grid is initialized 
  */
  init_grid(1<<maxlevel);
  const face vector muc[]={vis,vis};
  mu=muc;
  fprintf(ferr,"Re = %g and max res = %d X %d\n",1/vis,(int)pow(2,maxlevel),(int)pow(2,maxlevel));
  foreach()
    u.x[]=0.;
  boundary(all);
  DT=0.01;
}
/**
For $Re=500$ we output an animation of the vorticity and the grid resolution.  
*/
event movie(t+=0.1){
  if (m==2){
    fprintf(ferr,"%g\n",t);
    foreach()
      omg[]=((u.x[0,1]-u.x[0,-1]) - (u.y[1,0]-u.y[-1,0]))/(2*Delta);
    boundary({omg}); //For interpolation purposes
    static FILE * fp3 = popen ("ppm2gif > omg.gif", "w");
    output_ppm(omg,fp3, n=512,linear=true, min=-10., max=10.);
  }
}

/**
At $t=20L/U$ the simulation is stopped and some output statistics are calculated.
*/
event stop(t+=1.;t<25.){
  if (t>=20.){
    ens=0;
    vit=0;
    vir=0;
    e=0;
    foreach(reduction(+:ens) reduction(+:e)){
      omg[]=((u.x[0,1]-u.x[0,-1]) - (u.y[1,0]-u.y[-1,0]))/(2*Delta);
      ens+=sq(omg[])*sq(Delta);
      lev[]=level;
      U[]=sqrt(sq(u.x[])+sq(u.y[]));
      e+=sq(U[])*sq(Delta);
    }
    boundary({omg});
    foreach_boundary(top)
      vit+=Delta*omg[];
    foreach_boundary(right)
      vir+=Delta*omg[];
    foreach_boundary(left)
      vil+=Delta*omg[];
    foreach_boundary(bottom)
      vib+=Delta*omg[];
    fprintf(ferr,"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",i,t,e,ens,vit,vir,vir,vib);
    return 1;
  }
}

/**
## Results
We look at the evolution of the voricity field for $Re=500$ 

![Evolution of the vorticity field](lidfg/omg.gif)

*/