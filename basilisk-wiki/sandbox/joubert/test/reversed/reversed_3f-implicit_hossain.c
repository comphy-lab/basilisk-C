/**
# Three fluids time-reversed VOF advection in a vortex.

This is a modified version of reversed.c test which advects and stretches an
initially circular interface composed of two fluid fractions in a
non-divergent vortical flow in third fluid. The flow reverts in time and the interface
should come back to its original position. The difference between the initial
and final shapes is a measure of the errors accumulated during advection.

We will need the advection solver combined with the VOF advection scheme. */

#include "advection.h"
/**
  The volume fraction are stored in scalar field `f1`,`f2`,`f3` which are listed as
  an *interface* for the VOF solver. We do not advect any tracer with
  the default (diffusive) advection scheme of the advection solver. */

scalar f1[], f2[], f3[];
scalar * interfaces = {f1,f2,f3}, * tracers = NULL;

int MAXLEVEL;
double volref[3], sum;
FILE * fp2 = NULL;
FILE * fp3 = NULL;
FILE * fp1 = NULL;
FILE * fp = NULL;

#if Normalize
// We include regular Vof function.
#include "vof.h"
#else
//We include the three phase vof function designed by Hossain Chizari.
#ifndef R_VOFLIMIT
#define R_VOFLIMIT 1e-3
#endif
#include "../../chizari/threephase/vof-3p.h"
#endif

/**
  We center the unit box on the origin and set a maximum timestep of 0.1 */

int main() {
  origin (-0.5, -0.5);
  DT = .1;

  fp3 = fopen ("f3", "w");
  fp2 = fopen ("f2", "w");
  fp1 = fopen ("diff", "w");
  fp = fopen ("mass", "w");

  /**
    We then run the simulation for different levels of refinement. */

  for (MAXLEVEL = 5; MAXLEVEL <= 8; MAXLEVEL++) {
    init_grid (1 << MAXLEVEL);
    run();
    fprintf (stderr,"\n\n");
    fprintf (fp,"\n\n");

  }
  fclose (fp);
  fclose (fp1);
  fclose (fp2);
  fclose (fp3);
}

/**
  The initial interface is a circle of radius 0.2 centered on
  (-0.236338,-0.2) (for historical reasons). We use the levelset
  function `circle()` to define this interface. 

  The period of the stretching cycle is set to 15, which will lead to
  strong stretching. Milder conditions can be obtained by decreasing it. */

#define circle(x,y) (sq(0.2) - (sq(y + 0.2) + sq(x + .236338)))
#define T 15

/**
  We use the functions below to respectively define the bottom and upper halfs 
  of the circle and the outside using boolean operations.*/

double halfb (double x, double y){
  double lineb = (-0.236338-x);
  return min (circle(x,y),lineb);
}

double halfu (double x, double y){
  double lineu = (x+0.236338);
  return min (circle(x,y),lineu);
}

double outside (double x, double y){
  double blockb = (-y-0.0);
  double blocku = (y+0.0);
  return min (max(blockb,blocku), -circle(x,y));
}

/** 
  We compute the corresponding volume fraction field to the boolean
  geometry.*/

#if NORMALIZE
event init (i = 0) {
  fraction (f1, halfb(x,y));
  fraction (f2, halfu(x,y));
  fraction (f3, outside(x,y));

  volref[0] = statsf(f1).sum;
  volref[1] = statsf(f2).sum;
  volref[2] = statsf(f3).sum;

#else
  event init (i = 0) {
    fraction (f1, halfb(x,y));
    fraction (f2, halfu(x,y));

    foreach()
      f3[] = 1.0 - f1[] - f2[];

    volref[0] = statsf(f1).sum;
    volref[1] = statsf(f2).sum;
    volref[2] = statsf(f3).sum;
  }
#endif

  event properties (i++){
    sum=0.0;
    foreach ()
    {
      f3[] = 1.0 - f1[] - f2[];
      if (f1[] + f2[] > 1.0 - R_VOFLIMIT)
      {
	sum = f1[] + f2[];
	f3[] = fabs(0.0);
	f2[] /= sum;
	f1[] /= sum;
      }
      else if (f1[] + f2[] < R_VOFLIMIT)
      {
	f3[] = 1.0;
	f2[] = fabs(0.0);
	f1[] = fabs(0.0);
      }
      else
      {
	if (f2[] < 0.0 + R_VOFLIMIT)
	{
	  sum = f3[] + f1[];
	  f1[] /= sum;
	  f2[] = fabs(0.0);
	  f3[] /= sum;
	}
	else if (f2[] > 1.0 - R_VOFLIMIT)
	{
	  f1[] = fabs(0.0);
	  f2[] = 1.0;
	  f3[] = fabs(0.0);
	}
	if (f1[] < 0.0 + R_VOFLIMIT)
	{
	  sum = f2[] + f3[];
	  f3[] /= sum;
	  f2[] /= sum;
	  f1[] = fabs(0.0);
	}
	else if (f1[] > 1.0 - R_VOFLIMIT)
	{
	  f3[] = fabs(0.0);
	  f2[] = fabs(0.0);
	  f1[] = 1.0;
	}
      }
    }
    boundary({f1, f2, f3});
  }

#if NORMALIZE
  event fnormalize (i++){
    sum=0.0;
    foreach(){
      sum = (f1[]+f2[]+f3[]);
      f1[]= clamp((f1[]/sum),0.,1.);
      f2[]= clamp((f2[]/sum),0.,1.);
      f3[]= clamp((f3[]/sum),0.,1.);
    }
    boundary ({f1,f2,f3});
  }
#endif

  event velocity (i++) {

    /**
      This event defines the velocity field.  
      On trees we first adapt the grid so that the estimated error on
      the volume fractions is smaller than $5\times 10^{-3}$. We limit the
      resolution at `MAXLEVEL` and we refine the volume fraction fields
      `f1` and `f2`. */
#if TREE

    adapt_wavelet ({f1,f2,f3}, (double[]){5e-3,5e-3,5e-3}, MAXLEVEL, list = {f1,f2,f3});

#endif

    /**
      The velocity field is defined through a streamfunction $\psi$, defined
      on the vertices of the grid. */

    vertex scalar psi[];
    foreach_vertex()
      psi[] = - 1.5*sin(2.*pi*t/T)*sin((x + 0.5)*pi)*sin((y + 0.5)*pi)/pi;

    /** 
      We can then differentiate the streamfunction to get the velocity components.
      This guarantees that the velocity field is exactly non-divergent. */

    trash ({u});
    struct { double x, y; } f = {-1.,1.};
    foreach_face()
      u.x[] = f.x*(psi[0,1] - psi[])/Delta;
    boundary ((scalar *){u});

  }

  /**
    At the start and end of the simulation we check the sum, min and max values
    of the sum of the two volume fraction field. The sum must be constant to
    within machine precision and the total volume fraction should be bounded by
    zero and one. */

  event massfiles (t += T/15.) {
    scalar  mass[];
    foreach()
      mass[] =(f1[]+f2[]+f3[]);

    stats s = statsf (mass);
    fprintf (fp, "%g %g %g %g %g\n", t, s.sum, s.min, s.max,
	normf (mass).rms);
  }

  event logfile (i++)
  {
    fprintf (stderr, "%g %g %g %g %g\n", t,
	statsf(f1).sum - volref[0],
	statsf(f2).sum - volref[1],
	statsf(f3).sum - volref[2],
	normf (u.x).max);
  }

  /**
    To compute the error, we reinitialise field `e` at the end of the
    simulation with the initial circle shape and compute the difference with the
    final shape. We output the norms as functions of the maximum
    resolution `N`. */

  event field (t = T) {
    scalar e[], k[], j[];

    fraction (e, halfb(x,y));
    fraction (k, halfu(x,y));
    fraction (j, outside(x,y));

    foreach(){
      e[] -= f1[];
      k[] -= f2[];
      j[] -= f3[];
    }
    norm n = normf (e);
    norm m = normf (k);
    norm o = normf (j);
    fprintf (fp1, "%d %g %g %g %g %g %g %g %g %g\n", N, n.avg, n.rms, n.max, 
	m.avg, m.rms, m.max, o.avg, o.rms, o.max);
  }

  /**
    We also output the shape of the reconstructed interfaces at regular
    intervals (but only on the finest grid considered). */
  event shape (t += T/4.) {
    if (N == 256){
      output_facets (f1, stdout);
      output_facets (f2, fp2);
      output_facets (f3, fp3);
    }
  }

  /**
    If we are using adaptivity, we also output the levels of refinement at
    maximum stretching. */

#if TREE

  event levels (t = T/2) {
    if (N == 256) {
      scalar l[];
      foreach()
	l[] = level;
      output_ppm (l, file = "levels.png", n = 400, min = 0, max = 7);
    }
  }

#endif

#if 0

  event movie (i += 10)
  {
    scalar l[];
    foreach()
      l[] = level;
    output_ppm (l, 1 << MAXLEVEL, file = "level.ppm");
  }

#endif

  /**

##Results

~~~gnuplot Shapes of the interface for $t=0$, $t=T/4$, $t=T/2$, $t=3T/4$ and $t=T$ 
reset
set xlabel 'x'
set ylabel 'y'
set size ratio -1
plot 'out' w l t 'f1', 'f2' w l t 'f2',\
'f3' w l t 'f3'
~~~

~~~gnuplot Convergence of error between initial and final shape in function of cells number
reset
set xlabel 'Number of cells'
set ylabel 'Maximum error'
set yrange [0.01:10]
set logscale x
set logscale y
plot 'diff' index 0:0 u 1:3 w l t 'error(f1) rms', \
'diff' index 0:0 u 1:6 w l t 'error(f2) rms', \
'diff' index 0:0 u 1:9 w l t 'error(f3) rms'
~~~

~~~gnuplot Mass conservation (volume sum) in function of time
reset
set xlabel 'time'
set ylabel 'Mass(f1[]+f2[]+f3[]) sum'
plot 'mass' index 0:0  u 1:2 w l t 'f1+f2+f3 lvl5',\
'mass' index 1:1  u 1:2 w l t 'f1+f2+f3 lvl6',\
'mass' index 2:2  u 1:2 w l t 'f1+f2+f3 lvl7',\
'mass' index 3:3  u 1:2 w l t 'f1+f2+f3 lvl8'
~~~

~~~gnuplot Mass conservation (RMS norm) in function of time
reset
set xlabel 'time'
set ylabel 'Mass(f1[]+f2[]+f3[]) RMS norm'
plot 'mass' index 0:0  u 1:5 w l t 'f1+f2+f3 lvl5',\
'mass' index 1:1  u 1:5 w l t 'f1+f2+f3 lvl6',\
'mass' index 2:2  u 1:5 w l t 'f1+f2+f3 lvl7',\
'mass' index 3:3  u 1:5 w l t 'f1+f2+f3 lvl8'
~~~

~~~gnuplot Convergence of error (volume sum) on f1 in function of time
reset
set xlabel 'time'
set ylabel 'Error on f1'
plot 'log' index 0:0  u 1:2 w l t 'f1 lvl5',\
'log' index 1:1  u 1:2 w l t 'f1 lvl6',\
'log' index 2:2  u 1:2 w l t 'f1 lvl7',\
'log' index 3:3  u 1:2 w l t 'f1 lvl8'
~~~

~~~gnuplot Convergence of error (volume sum) for one level in function of time
reset
set xlabel 'time'
set ylabel 'Error on fluid fraction'
plot 'log' index 3:3  u 1:2 w l t 'f1 lvl8',\
'log' index 3:3  u 1:3 w l t 'f2 lvl8',\
'log' index 3:3  u 1:4 w l t 'f3 lvl8'
~~~
*/
