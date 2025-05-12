/**
# Direct numerical simulation of convective turulence
In this test case adaptive grids are tested for convective turbulence with a fixed bottom temperature (i.e. Buoyancy) and (inital) liniarly stratified "atmosphere". The direct numerical simulation (DNS) test case is inspired by [1], Here reference results are provided, the used code here is [MicroHH](http://microhh.org/).

## Loading modules and declaring parameters
- We specify that we would like to perform our computation on an oct-tree grid (i.e. 3D adaptive) 
- We specify that we want to solve the equations of fluid motion, specify a diffusive tracer (T) that will be our buoyancy field
- Set the maximum level of refinment (minimum is " MAXLEVEL - 4" ) 
- We initizalize some global variables that will prove to be usefull.

Note: We have a normalized linear stratification, bottom buoyancy and thereby lengthscale L. Therefore, the Reynoldsnumber is directly controlled by the viscosity of the fluid. 
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "profil/profile2.h"
#include "lambda2.h"
#define sqN 1
#define MAXLEVEL 9
#define TT pow(molvis,-0.333333333)
#define temp 45*TT

#define sq(x)  x * x 

char names[80];   
double ue = 1e-2; 
scalar T[];
scalar * tracers = {T};
mgstats mgT;
face vector av[];

double molvis = 1.225289152e-4;

/**
## Main
Set boundary conditions, define domain, initialize viscosity, intialize grid, link boussinesq accelertation term, and call "run()" loop for time integration.
*/
int main() 
{
  T[bottom]=dirichlet(1);
  T[top] = dirichlet(3);
  u.t[bottom]=dirichlet(0);
  L0 = 3.;
  X0 = Y0 = Z0 = 0;
  init_grid(1 << (MAXLEVEL-1));
  const face vector muc[] = {molvis,molvis,molvis};
  mu = muc;
  a=av;
  run();
}
/**
## Initialization 
- Set a small timestep and tolerance on the poisson solver. This is for the first few timesteps for a correct initisalization
- Initizalize the stratification
- Refine the grid and add a perturbation to buoyancy and velocity components to kickstart the growth of the instability. 
*/

event init(t=0)
{      
    DT=0.001;
  TOLERANCE=1e-5; 
  foreach()
    {
      T[]=(sqN*y);
    }
  refine (y < 0.05 && level < (MAXLEVEL),  all);

  foreach()
    {
      T[]+=(0.04*noise()*(exp(-y/0.02)));
      foreach_dimension()
	{
	  u.x[]=0.01*noise()*(exp(-y/0.02));
	}
    }
  boundary(all);
}
 
/** 
## Buoyancy acceleration & Spongelayer and diffusion of T[] field
Define the acceleration term and do diffusion of diffusive buoyancy scalar.
*/ 
event acceleration (i++)
{   
  coord Del = { 0 , 1, 0}; 
  foreach_face()
    {
      av.x[] = Del.x* (((T[]+T[-1])/2)-(((u.x[]+u.x[-1])/2)*((exp((y-3)/0.1)*(y>3)*(y<3))+(y>3)))); 
    }
  
  foreach()
    {
      T[]-=(T[]-y) * ((exp((y-3)/0.05)*(y>2)*(y<3))+(y>3));
  }
}
/** 
## Diffusion
Diffusion of buoyancy scalar.
*/ 
event Diffusion (i++)
{
  mgT = diffusion (T, dt, mu);
  boundary(all);
}
/**
## Log file
Boring log file
*/
event logfile(t+=1; t<=temp)
{
  
  fprintf(ferr,"%dD\tt=%g\ti=%d\tDT=%g\n" ,dimension,t,i,DT);
  
}
/**
## Adaptation 
Refine and coarsen the grid and set the timestep and tolerance to a more sensible value after some initial timesteps. Note however that, from this point on, the timestep is typically limited by the CFL condition, not by DT. 
*/ 
event adapt (i++)
{
  adapt_wavelet((scalar *){u,T},(double[]){ue,ue,ue,ue},MAXLEVEL,(MAXLEVEL-4));
  if (i==10)
  {
    fprintf(ferr,"Adapt timestep before DT = %g.",DT);
    DT=0.1;
    TOLERANCE = 1e-3; 
    fprintf(ferr,"now DT = %g\n",DT);
  }
}
/**
# Diagnostics
The code below is used to diagnose the evolution of the simulation and numerically obtained results with vertical profiles, time series and slices.

## Profiles
First we use the 'profile()' function included from profile2.h to write files containing vertical (y-direction) profiles of the vertical variance, the horizontal variance, the kinetic energy, the resolved dissipation and the ratio between an estimated viscous length scale and the grid box size. The profiels are evaluated for [0<y<2] from slices taken at the maximum resolution. 
*/

event profiler(t+=TT)
{
  scalar vervar[];
  scalar horvar[];
  scalar energy[];
  scalar dis[];
  scalar kolrat[];
  foreach()
    {
      vervar[] = sq(u.y[]);
      horvar[] = sq(u.x[])+sq(u.z[]); 
      energy[] = sq(u.x[])+sq(u.y[])+sq(u.z[]);
      dis[]  = molvis*
	(sq(((u.x[1] - u.x[-1])/(2*Delta))) +
	 sq(((u.x[0,1] - u.x[0,-1])/(2*Delta))) +
	 sq(((u.x[0,0,1] - u.x[0,0,-1])/(2*Delta)))+
	 sq(((u.y[1] - u.y[-1])/(2*Delta))) +
	 sq(((u.y[0,1] - u.y[0,-1])/(2*Delta))) +
	 sq(((u.y[0,0,1] - u.y[0,0,-1])/(2*Delta)))+
	 sq(((u.z[1] - u.z[-1])/(2*Delta))) +
	 sq(((u.z[0,1] - u.z[0,-1])/(2*Delta))) +
	 sq(((u.z[0,0,1] - u.z[0,0,-1])/(2*Delta))));
    }
  foreach()
    {
      kolrat[] = Delta*pow((dis[]/pow(molvis,3)),0.25);
    }
  boundary({T,energy,dis,horvar,vervar,kolrat});
  profile({T,energy,dis,horvar,vervar,kolrat},0, 2,(MAXLEVEL),t);
}

/**
## Time serie variables
Output variables for energy, Boundary layer height, number of grid cells, bottom boundary buoyancy flux (i.e. conductive), 
dissipation using different formulations, number of grid cells at certain levels, buoyancy flux (i.e. advective) and the ratio of an estimated viscous lengthscale and the grid box size. 
We integrate the domain as a discrete summation.   
*/
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
    hs += dv()*  (T[]-y);
  return hs;
}

static int n()
{
  int ne=0; 
  foreach(reduction(+:ne))
    ne+=1;
  return ne;
}

static double bf()
{
  double bufx = 0;
  foreach_boundary(bottom reduction(+:bufx))
    bufx+=sq(Delta)*molvis*(T[ghost] - T[])/Delta;
  return bufx;
}

static double diss()
{  double dis = 0, vol = 0;
  // Loop over iedere gridcell
  foreach(reduction(+:vol) reduction(+:dis))
    {
      vol += dv();
      
	// bereken voor iedere snelheidsrichting de afgeleiden in 3 richting richtingen, kwadrateer, tel op, en weeg met volume van de gridcell.
      
	  dis += dv()*(sq(((u.x[1] - u.x[-1])/(2*Delta))) +
		       sq(((u.x[0,1] - u.x[0,-1])/(2*Delta))) +
		       sq(((u.x[0,0,1] - u.x[0,0,-1])/(2*Delta)))+
		       sq(((u.y[1] - u.y[-1])/(2*Delta))) +
		       sq(((u.y[0,1] - u.y[0,-1])/(2*Delta))) +
		       sq(((u.y[0,0,1] - u.y[0,0,-1])/(2*Delta)))+
		       sq(((u.z[1] - u.z[-1])/(2*Delta))) +
		       sq(((u.z[0,1] - u.z[0,-1])/(2*Delta))) +
		       sq(((u.z[0,0,1] - u.z[0,0,-1])/(2*Delta))));
	
    }
  // Vermenigvildig met viscositeit. 
  dis *= molvis;
  return dis;
}

static double diss2()
{  double dis = 0, vol = 0;
  // Loop over iedere gridcell
  foreach(reduction(+:vol) reduction(+:dis))
    {
      vol += dv();
      
	// bereken voor iedere snelheidsrichting de afgeleiden in 3 richting richtingen, kwadrateer, tel op, en weeg met volume van de gridcell.
    
      dis += dv()*(sq(((u.x[1] - u.x[-1])/(2*Delta))) +
		   sq(((u.x[0,1] - u.x[0,-1])/(2*Delta))) +
		   sq(((u.x[0,0,1] - u.x[0,0,-1])/(2*Delta)))+
		   sq(((u.y[1] - u.y[])/(Delta))) +
		   sq(((u.y[0,1] - u.y[])/(Delta))) +
		   sq(((u.y[0,0,1] - u.y[])/(Delta)))+
		   sq(((u.z[1] - u.z[-1])/(2*Delta))) +
		   sq(((u.z[0,1] - u.z[0,-1])/(2*Delta))) +
		   sq(((u.z[0,0,1] - u.z[0,0,01])/(2*Delta))));
	
  
    }
  // Vermenigvildig met viscositeit. 
  dis *= molvis;
  return dis;
}

static double diss3()
{  double dis = 0, vol = 0;
  // Loop over iedere gridcell
  foreach(reduction(+:vol) reduction(+:dis))
    {
      vol += dv();
      
	// bereken voor iedere snelheidsrichting de afgeleiden in 3 richting richtingen, kwadrateer, tel op, en weeg met volume van de gridcell. 
      dis += dv()*(sq(((uf.x[] - uf.x[-1])/(Delta))) +
		   sq(((uf.x[] - uf.x[0,-1])/(Delta))) +
		   sq(((uf.x[] - uf.x[0,0,-1])/(Delta)))+
		   sq(((uf.y[] - uf.y[-1])/(Delta))) +
		   sq(((uf.y[] - uf.y[0,-1])/(Delta))) +
		   sq(((uf.y[] - uf.y[0,0,-1])/(Delta)))+
		   sq(((uf.z[] - uf.z[-1])/(Delta))) +
		   sq(((uf.z[] - uf.z[0,-1])/(Delta))) +
		   sq(((uf.z[] - uf.z[0,0,-1])/(Delta))));
  
    }
  // Vermenigvildig met viscositeit. 
  dis *= molvis;
  return dis;
}

static int n1()
{
  int nee = 0;
  foreach(reduction(+:nee))
    if (pid()==1)
      nee++;
  return nee;
}

static double bufx()
{
  double buff = 0; 
  double bg=0;
  for (double yp=L0/(pow(2,MAXLEVEL+1));yp<2;yp+=L0/(pow(2,MAXLEVEL)))
    {
      
      int i=0;
      
      bg = 0;
      double ** field = matrix_new (pow(2,MAXLEVEL),pow(2,MAXLEVEL),sizeof(double));
      double ** fiel = matrix_new (pow(2,MAXLEVEL),pow(2,MAXLEVEL),sizeof(double)); 
      for (double xp=L0/(pow(2,MAXLEVEL+1));xp < L0;xp+=L0/(pow(2,MAXLEVEL)))
	{
	  int k=0;
	   for (double zp=L0/(pow(2,MAXLEVEL+1));zp < L0;zp+=L0/(pow(2,MAXLEVEL)))
	     {
	       Point point = locate(xp,yp,zp); 
	       field[i][k] = point.level >= 0 ? interpolate(T,xp,yp,zp) : nodata;
	       fiel[i][k] = point.level >=0 ? interpolate(u.y,xp,yp,zp) : nodata;
	       //fprintf(ferr,"%g\t%d\t%d\n",field[i][k],i,k);
	      k++;
	     }
	   i++;
	}
      if (pid() == 0)
	{
	  // master
	  @ if _MPI
	    {
	      MPI_Reduce (MPI_IN_PLACE, field[0], pow(2,(MAXLEVEL*2)), MPI_DOUBLE, MPI_MIN, 0,
			MPI_COMM_WORLD);
	      MPI_Reduce (MPI_IN_PLACE, fiel[0], pow(2,(MAXLEVEL*2)), MPI_DOUBLE, MPI_MIN, 0,
		      MPI_COMM_WORLD);
	    }
	  @endif
	    for(int i = 0; i<pow(2,MAXLEVEL);i++)
	      {
		for (int k = 0; k<pow(2,MAXLEVEL);k++)
		  {
		    bg+=field[i][k];
		  }
	      }
	  bg = bg/pow(2,(MAXLEVEL*2)); 
	  
	  for (int i = 0; i<pow(2,MAXLEVEL);i++)
	    {
	      for (int k = 0;k<pow(2,MAXLEVEL);k++)
		{
		  buff+= (fiel[i][k]*(field[i][k]-bg))*pow(L0/pow(2,MAXLEVEL),3); 
		}
	    }
	  
	}
      @if _MPI
      else // slave
	{
	  MPI_Reduce (field[0], NULL, pow(2,MAXLEVEL*2), MPI_DOUBLE, MPI_MIN, 0,
		    MPI_COMM_WORLD);
	  MPI_Reduce (fiel[0], NULL, pow(2,MAXLEVEL*2), MPI_DOUBLE, MPI_MIN, 0,
		  MPI_COMM_WORLD);
	}
      @endif
	
	matrix_free (field);
        matrix_free (fiel); 
    }
  return buff;
}

static double kolrat()
{
  scalar dis[];
  scalar kolrt[];
  foreach()
    {
      dis[]  = (molvis)*
	(sq(((u.x[1] - u.x[-1])/(2*Delta))) +
	 sq(((u.x[0,1] - u.x[0,-1])/(2*Delta))) +
	 sq(((u.x[0,0,1] - u.x[0,0,-1])/(2*Delta)))+
	 sq(((u.y[1] - u.y[-1])/(2*Delta))) +
	 sq(((u.y[0,1] - u.y[0,-1])/(2*Delta))) +
	 sq(((u.y[0,0,1] - u.y[0,0,-1])/(2*Delta)))+
	 sq(((u.z[1] - u.z[-1])/(2*Delta))) +
	 sq(((u.z[0,1] - u.z[0,-1])/(2*Delta))) +
	 sq(((u.z[0,0,1] - u.z[0,0,-1])/(2*Delta))));
    }
  foreach()
    {
      kolrt[] = Delta*pow((dis[]/pow(molvis,3)),0.25);
    }
  stats kol = statsf(kolrt);
  return kol.max;
}
char namm[50];
event timeseries(t+=(TT/20))
{
  char nja[10]="./data"; 
  sprintf(namm,"%s/timeseries.dat",nja);
  static FILE * fp = fopen ("./data/timeseries", "w"); 
  if (t==0)
    fprintf(fp,"t\ttwall\tspeed\tdt\ti\tn\te\thst\tsbufx\tdissipation\tDissipation\tDISSIPATION\tkoldelrat\tn1\tbufx\n");
  fprintf(fp,"%g\t%g\t%g\t%g\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%g\n",t,perf.t,perf.speed,dt,i,n(),energy(),hst(),bf(),diss(),diss2(),diss3(),kolrat(),n1(),bufx()); 
  fflush(fp); 
   
}
/**
## Slices
Files Containing vertical and Horizontal slices of various interesting fields are writen by the event below. 
*/

event slices(t+=(TT))
{

  // slice of temperature- and gradient field
  scalar grT[];
  scalar piz[];
  scalar bufx[];
  scalar dis[];
  
  foreach()
    {
      grT[]= pow(sq((T[1]-T[-1])/ (2*Delta))+sq((T[0,1]-T[0,-1])/ (2*Delta))+sq((T[0,0,1]-T[0,0,-1])/ (2*Delta)),0.5);
      piz[] = tid();
      dis[]  = (molvis)*
	(sq(((u.x[1] - u.x[-1])/(2*Delta))) +
	 sq(((u.x[0,1] - u.x[0,-1])/(2*Delta))) +
	 sq(((u.x[0,0,1] - u.x[0,0,-1])/(2*Delta)))+
	 sq(((u.y[1] - u.y[-1])/(2*Delta))) +
	 sq(((u.y[0,1] - u.y[0,-1])/(2*Delta))) +
	 sq(((u.y[0,0,1] - u.y[0,0,-1])/(2*Delta)))+
	 sq(((u.z[1] - u.z[-1])/(2*Delta))) +
	 sq(((u.z[0,1] - u.z[0,-1])/(2*Delta))) +
	 sq(((u.z[0,0,1] - u.z[0,0,-1])/(2*Delta))));
    }  
  scalar * list = {T,grT,piz,dis};
  sprintf(names,"./data/verslicet=%g.dat",t);
  FILE *fpver =fopen (names,"w"); 
  int nn = (pow(2,MAXLEVEL));
  int len = list_len(list);
  double ** field = matrix_new (nn, nn, len*sizeof(double));
  double zp = L0/2.01;
  
  double stp = L0/nn;
  for (int i = 0; i < nn; i++)
    {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) 
	{
	  double yp = stp*j + Y0 + stp/2.;
	  Point point = locate (xp, yp,zp);
	  int k = 0;
	  for (scalar s in list)
	    field[i][len*j + k++] = point.level >= 0 ? s[] : nodata;
	  //field[i][len*j + k++] = interpolate (s, xp, yp,zp);
	}
    }
  if (pid() == 0)
    { // master
      @if _MPI
	MPI_Reduce (MPI_IN_PLACE, field[0], len*nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		    MPI_COMM_WORLD);
      @endif
	int k = 0;
      for (scalar s in list)
	{
	  for (int i = 0; i < nn; i++) {
	    for (int j = 0; j < nn; j++) {
	      fprintf (fpver, "%g\t", field[i][len*j + k]);
	    }
	    fputc ('\n', fpver);
	  }
	  fflush (fpver);
	  k++;
	}
    }
  @if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, len*nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
  @endif
    matrix_free (field);


  
  scalar * lists = {T,u.y};
  sprintf(names,"./data/horslicet=%g.dat",t);
  FILE *fphor =fopen (names,"w"); 
  nn = (pow(2,MAXLEVEL));
  len = list_len(lists);
  field = matrix_new (nn, nn, len*sizeof(double));
  double yp = Y0+(1/3);
  
  for (int i = 0; i < nn; i++)
    {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) 
	{
	  double zp = stp*j + Y0 + stp/2.;
	  Point point = locate (xp, yp,zp);
	  int k = 0;
	  for (scalar s in lists)
	    // field[i][len*j + k++] = interpolate (s, xp, yp,zp);
	  field[i][len*j + k++] = point.level >= 0 ? s[] : nodata;
	}
    }
  if (pid() == 0)
    { // master
      @if _MPI
	MPI_Reduce (MPI_IN_PLACE, field[0], len*nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		    MPI_COMM_WORLD);
      @endif
	int k = 0;
      for (scalar s in lists)
	{
	  for (int i = 0; i < nn; i++) {
	    for (int j = 0; j < nn; j++) {
	      fprintf (fphor, "%g\t", field[i][len*j + k]);
	    }
	    fputc ('\n', fphor);
	  }
	  fflush (fphor);
	  k++;
	}
    }
  @if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, len*nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
  @endif
    matrix_free (field);
}




/**
# Results 
A version of the code above was run on the SURF-sara HPC-cloud computer. A visualization of the results looks something like this: 
[Youtube link](https://www.youtube.com/watch?v=K1zTd5eoFQ4) 

A more quantitative validation with reference results is the comparison of the energy in the flow field as a function of time: 

![Validation results: adaptive solver (Basilisk) and MicroHH by reference[1]](/sandbox/Antoonvh/chielzor.jpg)

# Reference
[1] van Heerwaarden, Chiel C., and Juan Pedro Mellado. "Growth and decay of a convective boundary layer over a surface with a constant temperature." Journal of the Atmospheric Sciences 2016 (2016).
 */
