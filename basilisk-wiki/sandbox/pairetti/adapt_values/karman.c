/**
# Adapt_values function tutorial

This simple case is based on the [von Karman vortex street example](/src/examples/karman.c)
using the solid representation from the [moving cylinder case](/sandbox/popinet/moving_cylinder.c).

The aim of this tutorial is showing the [adapt_values function](http://basilisk.fr/sandbox/pairetti/adapt_values/adapt_values.h),
and how to use it in order to produce a *variable emax* simulation.

We first define the problem parameters and refinement settings.

*/

int MAXLEVEL = 9;
double emax = 0.05;
#define radius 0.0625
#define RE 160.
double maxRunTime = 20.;

/**
Flow viscosity is computed from radius and Reynolds. Density and inlet velocity
values are equal to 1.
*/
double muVal = 2*radius/RE;


/**
We are going to solve incompressible Navier-Stokes on collocated grid,
using a solid volume fraction to represent the cylinder.
*/
#include "navier-stokes/centered.h"
#include "fractions.h"
face vector muv[];

/**
The main function accepts MAXLEVEL and emax redefinition.
*/

int main(int argc, char * argv[])
{
  if (argc > 1)
    MAXLEVEL = atoi (argv[1]);
  if (argc > 2)
    emax = atof (argv[2]);
    
/**
The domain is a full square of L = 64 radius. The cylinder is placed at 10
radius from the inlet. The initial mesh is set 4 levels coarser than MAXLEVEL.
*/
  L0 = 128.*radius;
  origin (-0.625, -L0/2.);
  N = pow(2.,MAXLEVEL-4);
  
  mu = muv;
  run();
}

/**
The fluid is injected on the left boundary with a unit velocity.
An outflow condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);


event init (t = 0) {

  refine (x > -0.25 && x < 0.25 &&
          y > -0.25 && y < 0.25 && level < MAXLEVEL + 1);

}

/**
Following the moving cylinder tutorial, the cylinder volume fraction
is saved in the cylinder field (which is global for post-processing purposes).
*/

scalar cylinder[];

event moving_cylinder (i++) {
  coord vc = {0.,0.};
  fraction (cylinder, - (sq(x - vc.x*t) + sq(y - vc.y*t) - sq(radius)));

  /**
  We then use this (solid) volume fraction field to impose the
  corresponding velocity field inside the solid (as a volume-weighted 
  average). */
  
  foreach()
    foreach_dimension()
      u.x[] = cylinder[]*vc.x + (1. - cylinder[])*u.x[];
      
  boundary ((scalar *){u});
  
}

/**
The adapt_values function applies refinement to cells with field values
above *emax*. To replicate the adapt_wavelet behaviour, we can pre-compute
the wavelet field of the velocity. Still, we can apply any filter to
reduce the *err* field value.

In this example, the weight is a ramped value, and the *errMeasure*
taken is equal to max(wav(u))*weight.

The errMeasure field is also defined out the event for post-processing matters.

*/
#include "./adapt-values.h"
#define max(a,b) ((a) > (b) ? (a) : (b))
scalar errMeasure[];

event adapt (i++)
{
  double x1 = 0.625, x2 = 2., slope, orig;
  x1 = 10.;  x2 = 12.;    // uncomment this line to remove weighting
  slope = 1/(x1 - x2);  orig = -x2*slope;
  
  scalar wavx[], wavy[], wavf[];
  wavelet(u.x, wavx);         wavelet(u.y, wavy);
  wavelet(cylinder, wavf);
  
  foreach()
  {
    double weight;
    if(x < x1)
      weight = 1;
    else
      weight = orig + slope*x > 0 ? orig + slope*x : 0;

    errMeasure[] = fabsf(weight*max(wavx[], wavy[]));
  }
     
  boundary({errMeasure});
     
  adapt_values ({wavf,errMeasure}, (double[]){1e-4, emax}, MAXLEVEL, 2);
}


/**
We output some convergence statistics... */

event logfile (i++; t <= maxRunTime)
{
  if(i==0)
  {
      double minDelta = L0/pow(2.,MAXLEVEL);
      double ppDiam = 2*radius/minDelta;
      fprintf (stderr, "# CASE params: Re = %g, mu = %g, rad = %g, ppDiam = %g \n", 
                        RE, muVal, radius, ppDiam);
      fprintf (stderr, "# CASE settings: MAXLEVEL = %d, emax = %g\n",
              MAXLEVEL, emax);

      fprintf (stderr, "# 1:i 2:t 3:mgp.i 4:mgu.i 5:grid->tn 6:perf.t 7:perf.speed\n");
  }
      
  fprintf (stderr, "%d %g %d %d %ld %g %g\n", 
                    i, t, mgp.i, mgu.i, grid->tn, perf.t, perf.speed);
}

/**
... and movies showing vorticity, refinement and errMeasure. 
Comparing, refinement and errMeasure you'll see that refinement
is limited to the zones where weight function is non zero.

![Animation of the vorticity field.](karman/vort.mp4)
![Animation of the refinement.](karman/level.mp4)
![Animation of the error function indicator.](karman/err.mp4)
*/

event movie (t += 0.1; t <= maxRunTime)
{

  scalar l[], m[], omega[];

  foreach()
  {
    l[] = level;
    m[] = 0.5 - cylinder[];
    omega[] = (u.y[1,0] - u.y[-1,0] + u.x[0,-1] - u.x[0,1])/(2.*Delta);
  }

  boundary ({omega});

  static FILE * fp = popen ("ppm2mp4 vort.mp4", "w");
  output_ppm (omega, fp, n = 1024, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m,
        min =  -50, max = 50 ,linear = false);

  static FILE * fp2 = popen ("ppm2mp4 level.mp4", "w");
  output_ppm (l, fp2, n = 1024, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m,
        min =  2, max = MAXLEVEL ,linear = false);
        
  static FILE * fp3 = popen ("ppm2mp4 err.mp4", "w");
  output_ppm (errMeasure, fp3, n = 1024, box = {{-0.5,-0.5},{7.5,0.5}},
        min =  0., max = emax ,linear = false);

}