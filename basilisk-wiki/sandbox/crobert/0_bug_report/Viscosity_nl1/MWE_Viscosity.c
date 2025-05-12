/** 
# Viscous diffusion
*/

#include "grid/multigrid1D.h"
//#include "layered/hydro.h"
#include "./hydro.h" //Uncomment to put the patch


/**

We test the ability of the [multilayer solver](/src/layered/hydro.h)
to resolve vertical viscosity in a flat parralel flow.

The interface is initialised with a uniform height, a unique cell in 
the horizontal direction and periodic boundary conditions.

Viscosity and total height are set to 1. Two cases are checked : 
* The Poiseuille flow with a viscous flow driven by a pressure gradient in the film
* The Couette flow with a viscous flow driven by a surfacic stress only (dut = 1).
*/
int ne = 0;
int surface = 0;
double a = 0.;
double emaxtrans = 0;

int main()
{
  linearised = true;
  N = 1;
  periodic(right);
  nu = 1.;
    
/** Poiseuille flow
*/  
  surface = 0;
  dut = zeroc;
  a = 1.;
  for (nl = 1; nl <= 8; nl += 1){
    ne = 0;
    run();
  }

/** Couette flow
*/
  surface = 1;
  dut = unity;
  a = 0.;
  for (nl = 1; nl <= 8; nl += 1){
    ne = 0;
    run();
  }
}

event init (i = 0)
{
  emaxtrans = 0;
  foreach()
    foreach_layer()
      h[] = 1./nl;
}

event pressure (i++)
{
  foreach_face()
    foreach_layer(){
   	  hu.x[] -= hf.x[]*dt*a;
	    ha.x[] -= hf.x[]*a;
    }
}

event timerfunction (t += 1e-5; t = 0; t <= 2.)
{}

#include "./ViscousDiffusion.h"

event logfile (t += 0.01; t = 0; t <= 2.)
{
  char name[80];
  sprintf (name, "velocity-nl%d-S%d", nl, surface);
  static FILE * fph = fopen (name, "w"); 
  fprintf(fph, "\n%g  ", t);
  double Q = 0;
  foreach()
    foreach_layer() {
      fprintf(fph, "%g ", u.x[]);
      Q += u.x[]*h[];
  }
  fprintf(fph, "   %g  %g  ", Q, Q - diffusion[ne][1]);
  emaxtrans = max (emaxtrans, Q - diffusion[ne][1]);
  ne++;
}

event error (t = end) {
  char name[80];
  sprintf (name, "final-nl%d-S%d", nl, surface);
  FILE * fpf = fopen (name, "w");
  foreach(){
    double z = 0.;
    foreach_layer(){
      z += h[]/2.;        
      fprintf(fpf, "\n %g  %g ", z, u.x[]);
      z += h[]/2.;  
    }      
  }
  

  if (surface == 1){
    static FILE * fp0 = fopen ("error-S0", "w");  
    fprintf(fp0, "%d %g \n", nl, emaxtrans);
  }
  else{
    static FILE * fp1 = fopen ("error-S1", "w");  
    fprintf(fp1, "%d %g \n", nl, emaxtrans);
  }
  
}


/**
# Results

#For the Poiseuille flow
~~~gnuplot Flux over time
set style line 1 pt 7 ps 0.7

set terminal @SVG enhanced size 640,640
plot 'velocity-nl1-S0' u 1:3 w l, \
     'velocity-nl2-S0' u 1:4 w l, \
     'velocity-nl3-S0' u 1:5 w l, \
     'velocity-nl4-S0' u 1:6 w l, \
     'velocity-nl5-S0' u 1:7 w l, \
     'velocity-nl6-S0' u 1:8 w l, \
     'velocity-nl7-S0' u 1:9 w l, \
     'velocity-nl8-S0' u 1:10 w l
~~~

~~~gnuplot Velocity over z at the end
plot 'final-nl1-S0' u 2:1 w p ps 3  ,\
     'final-nl2-S0' u 2:1 w l, \
     'final-nl3-S0' u 2:1 w l, \
     'final-nl4-S0' u 2:1 w l, \
     'final-nl5-S0' u 2:1 w l, \
     'final-nl6-S0' u 2:1 w l, \
     'final-nl7-S0' u 2:1 w l, \
     'final-nl8-S0' u 2:1 w l
~~~

#For the Couette flow
~~~gnuplot Flux over time
plot 'velocity-nl1-S1' u 1:3 w l, \
     'velocity-nl2-S1' u 1:4 w l, \
     'velocity-nl3-S1' u 1:5 w l, \
     'velocity-nl4-S1' u 1:6 w l, \
     'velocity-nl5-S1' u 1:7 w l, \
     'velocity-nl6-S1' u 1:8 w l, \
     'velocity-nl7-S1' u 1:9 w l, \
     'velocity-nl8-S1' u 1:10 w l
~~~

~~~gnuplot Error on the flux over time
plot 'velocity-nl1-S1' u 1:4 w l, \
     'velocity-nl2-S1' u 1:5 w l, \
     'velocity-nl3-S1' u 1:6 w l, \
     'velocity-nl4-S1' u 1:7 w l, \
     'velocity-nl5-S1' u 1:8 w l, \
     'velocity-nl6-S1' u 1:9 w l, \
     'velocity-nl7-S1' u 1:10 w l, \
     'velocity-nl8-S1' u 1:11 w l
~~~

~~~gnuplot Velocity over z at the end
plot 'final-nl1-S1' u 2:1 w p ps 3  ,\
     'final-nl2-S1' u 2:1 w l, \
     'final-nl3-S1' u 2:1 w l, \
     'final-nl4-S1' u 2:1 w l, \
     'final-nl5-S1' u 2:1 w l, \
     'final-nl6-S1' u 2:1 w l, \
     'final-nl7-S1' u 2:1 w l, \
     'final-nl8-S1' u 2:1 w l
~~~
*/