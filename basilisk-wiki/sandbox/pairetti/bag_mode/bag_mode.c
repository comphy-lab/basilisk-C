/**
# Droplet blow dy stream, bag mode secondary atomisation

We wish to study the behaviour of a single drop in free fall
affected dy a uniform flow of air, as in the [Opfer 2014](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.013023).

We use the centered Navier--Stokes solver. There are two phases, air and the liquid. The grid adaptation has been
modyfied to support region-depending maximum refinement (MAXLEVEL).*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#include "adapt_wavelet_limited.h"
bool limitedAdaptation = 0;

// bounding box values coordinates must be global
double xmin, xmax, ymin, ymax;

//Bounding box function (to define region of refinement)
void approx_bbox (scalar c){
  //gets the square where mixed cells are contained
  // saves it with xmin,ymin,xmax,ymax
  assert (dimension == 2);
  xmin = L0;
  xmax = -L0;
  ymin = L0;
  ymax = -L0;
  
  foreach()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      if(xmin > x)xmin = x;
      if(ymin > y) ymin = y;
      if(xmax < x) xmax = x;
      if(ymax < y) ymax = y;
    }
}

/**
Physical properties are defined based on [Opfer 2014](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.013023). */

/**Liquid properties*/
#define RHOL 785.
#define MUL 16.1e-3
#define SIGMA 48.3e-3

/**Air properties*/
#define RHOG 1.
#define MUG 1.81e-5

/**Drop diameter and stream velocity*/
#define DIAMETER 2.7e-3 
#define USTREAM 15.

#define MAXTIME 30e-3


/**
Minimum refinement level will be used directly, MAXLEVEL will only apply
on a certain region of the domain. */

int minlevel = 6;   // init grid of 32 points
int MAXLEVEL = 10;  //12

double maxruntime = HUGE;

/**
Boundary conditions are set only at left (inlet) and right (outlet)
patches. Inlet is supposed to be uniform */

u.n[left]  = dirichlet(USTREAM);
u.t[left]  = dirichlet(0);
u.n[right] = neumann(0);

p[left]    = neumann(0);
p[right]   = dirichlet(0);

/**
The main function can take the maximum refinement level as argument.
The main program will run ten steps with both versions of AMR.*/

int main (int argc, char * argv[]) {
  if (argc > 1)
    MAXLEVEL = atoi (argv[1]);
  
  /**
  We set the domain geometry and initial refinement.
  The drop will be centered at the origin, slightly
  moved to the left side*/
  L0 = 5*DIAMETER;
  size (L0);
  origin(-L0/4,-L0/2);
  init_grid (pow(2.0,minlevel));

  /**
  Time-step is limited to avoid stability issues related to surface tension
  [Galusinski & Vigneaux, 2008](http://www.sciencedirect.com/science/article/pii/S0021999108001253)*/
  double minDelta = L0/pow(2.0,MAXLEVEL);
  double sigmaTSLimit = sqrt(min(RHOL, RHOG)*pow(minDelta,3.0)/(SIGMA));
  double muTSLimit = min(MUL, MUG)*minDelta/SIGMA;
  DT = min(sigmaTSLimit,muTSLimit);

  // CFL number
  CFL = 0.4;

  /**
  Physical parameters are set by the macros at the beggining.*/
  
  f.sigma = SIGMA;
  rho1 = RHOL;
  rho2 = RHOG;

  mu1 = MUL;
  mu2 = MUG;

  /**
  Tolerance of the Poisson problem should never be reduced,
  this could lead to non divergence free velocity fields
  and mass conservation errors. */
  
  TOLERANCE = 1e-4; 

  run();

  limitedAdaptation = 1;
  run();
 
}

/**
Initial condition is given by drop position only, there is no flow at the begining.
First step should give a velocity field close to potential flow solution around the drop,
plus small velocities inside the droplet.*/

event init (t = 0) {
  if (!restore (file = "dump")){
    refine (sq(x) + sq(y) - sq(0.6*DIAMETER) < 0 && sq(x) + sq(y) - sq(0.4*DIAMETER) > 0 && level < MAXLEVEL);
    fraction (f, sq(0.5*DIAMETER) - (sq(x) + sq(y)));
  }
}

/**
Only ten steps are performed, just to check mesh adaptation.*/

event end (i = 10) {
  printf ("i = %d t = %g\n", i, t);
}

/**
## Region limited mesh adaptation

Any function with input(x,y,z) returning an integer can be used to
define MAXLEVEL locally for certain regions. Next, a simple example
using a circle (only for 2D case) slightly bigger than the drop is used.
Inside that region, AMR can refine until MAXLEVEL, outside
only MAXLEVEL-2 cells are allowed. Then, the alternative function
adapt_wavelet_limited is employed.*/

int refRegion(double x,double y, double z){
    int lev;
    if( sq(x)+sq(y) < sq(DIAMETER*0.7) )
      lev = MAXLEVEL;
    else
      lev = MAXLEVEL-4;

    return lev;
}

event adapt (i++) {
  double uemax = 5e-2;
  /**
  Choose between region limited adaptation or standard adaptation:*/
  if(limitedAdaptation)
    adapt_wavelet_limited ({f,u}, (double[]){1e-3,uemax,uemax,uemax}, refRegion, minlevel);
  else
    adapt_wavelet ({f,u}, (double[]){1e-3,uemax,uemax,uemax}, MAXLEVEL, minlevel);
}

/**
Every ten timesteps, we output the time, timestep, multigrid iterations, CPU time, speed and amount of cells. */

event logfile (i ++,first) {
  if (i == 0){
    fprintf (ferr,
	     "t dt mgp.i mgpf.i mgu.i perf.t perf.speed grid->tn\n");
  }

  fprintf (ferr,
	   "%g %g %d %d %d"
	   "%.2e %.2e %ld \n", 
	   t, dt, mgp.i, mgpf.i, mgu.i,
	   perf.t, perf.speed, grid->tn);
  fflush (ferr);
}

/**
First ten steps are saved in gfsview to compare initial flow field and mesh. */

event snapshot (i=0;i++; i <= 10)
{
   
  scalar omegaz[];
  foreach()
    omegaz[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
  boundary ({omegaz});
  
  char name[80];
  if(limitedAdaptation)
    sprintf (name, "snapshot-limited-%i_ts.gfs", i);
  else
    sprintf (name, "snapshot-standard-%i_ts.gfs", i);
  output_gfs (file = name, t = t, list = {f,u,p, omegaz});
}


/**
##Mesh evolution

From the log of this example, a graph as the following can be obtained:

[Mesh size during runtime](/bagMode_meshInTime.pdf)

*/