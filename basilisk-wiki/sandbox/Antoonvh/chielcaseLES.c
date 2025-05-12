/**
## Introducion 
This case is an example how to apply the Eddyviscosity function Vreman.h. that is located under "basilisk/src/SGS/ on my machine"
 The example case discrribes free convection into a liniarly stratified fluid with a fixed  bottom temperature

## step 1
Frist we include the octreegrid (to tell the simuliation it is 3D), eddyviscosity model, Navierstokes solver, and tracer and diffusion for buoyancy field
Next we define the stratification strength, maximum resolution level, time of simulation and the square operator
Finally, we declare some fields, values and strings that will be usefull during the simulation
 */
#include "grid/octree.h"
#include "SGS/Vreman.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

#define sqN 1
#define MAXLEVEL 10
#define temp 1000
#define sq(x)  x * x 

char names[80];   
scalar T[];
scalar * tracers = {T};
mgstats mgT;
face vector av[];
double molvis = 2.07841e-5;
face vector muv[];


/**
##Main 
set boundary conditions, domain size, initial grid, call the acceleration vector (for buoyancy and spongelayer) and the main timeloop */ 
int main() 
{
  T[bottom]=dirichlet(1);
  T[top] = dirichlet(3);
  u.t[bottom]=dirichlet(0);
  
  L0 = 6.;
  X0 = Y0 = 0;
  init_grid(1 << (MAXLEVEL-1));
  a=av;
  run();
}
/**
 ## initialize
call viscosity vector (to be determined later) initialize stratification and slice 50% of the top of the domain. 
*/
event init(t=0)
{
  mu=muv;
  mask(y>3 ? top : none);
  foreach()
    {
      T[]=(sqN*y);
    }
  boundary(all);
  
}
/** 
perturbate the stratification in order to trigger an instability to develop. This is done afther a few timesteps such that the grid at the bottom is refined at MAXLEVEL (a jump in temperature is located here initially).
 */ 
event perturb (i=4)
{
  foreach()
    {
      T[]+=(0.04*noise()*(exp(-y/0.02)));
    }
  boundary(all);
}

/*
Calculate Eddy viscosity (centered) scalar field ,turn it into a (face) vector field and set DT according to CFL condition for diffusion
*/

event SGS (i++)
{
  scalar muB[];
  eddyviscosity(0.17,u,muB,molvis);
  boundary({muB});
  foreach_face()
    {
      muv.x[]=(muB[]+muB[-1])/2;
    }
  boundary((scalar *){muv});

  scalar dft[];
  foreach()
    {	
      dft[] = 1*sq(Delta)/(muB[]);
    }	
  stats s = statsf(dft);
  DT= s.min;
}
/* apply boussinesq gravity and  sponge layer
 */
event acceleration (i++)
{   
  coord Del = { 0 , 1, 0}; 
  foreach_face()
    {
      av.x[] = Del.x* (((T[]+T[-1])/2)+((u.x[]*exp((y-3)/0.1)*(y>1.6)))); 
    }
  
  foreach()
    {
      T[]-=(T[]-y) * (exp((y-3)/0.05)*(y>2));
    }
  
}


event Diffusion (i++)
{
  mgT = diffusion (T, dt, mu);
  boundary(all);
}

event logfile(t+=1; t<=temp)
{
  fprintf(ferr,"%dD\tt=%g\ti=%d\tDT=%g\n" ,dimension,t,i,DT);
}
/* Adapt the grid. Currently i am looking into how the discretization error controls the gridsize. In a bad case scenario the simulation does refine toomuch and basically becomes a DNS.  
 */  
event adapt (i++)
{
  adapt_wavelet((scalar *){u,T},(double[]){4e-2,4e-2,4e-2},MAXLEVEL,6);
}


static double energy() 
{
  double e=0; 
  foreach(reduction(+:e))
    e += dv()* ( sq(u.x[]) + sq(u.y[]) + sq(u.z[])) ;
  return e;
}

static double hst() 
{
  double hs=0; 
  foreach(reduction(+:hs))
    hs += dv() *  (T[]-y);
  return hs;
}

static int n()
{
   int ne=0; 
   foreach(reduction(+:ne))
     ne+=1;
   return ne;
}

/*
  export flowfield quantities such as energy, boundarylayer height proxy (hst) and gridcell count  
*/

event timeseries(t+=1)
{
  if (t==0)
    {
      sprintf(names,"./datach/timeseriesLES.dat");
      FILE * fpp = fopen (names, "w");
      fprintf(fpp,"t\t\ti\t\t#n\t\te\t\thst\n");
      fclose(fpp);
    }
  FILE * fpp = fopen (names, "a");
  fprintf(fpp,"%g\t\t%d\t\t%d\t\t%g\t\t%g\n",t,i,n(),energy(),hst()); 
  fclose(fpp); 
}

/* 
Validation with DNS results is in progress...
 */ 
