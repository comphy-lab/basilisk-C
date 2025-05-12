/**
# Diffusion equation and Self Similar Solution(ML solver)
We want to solve the dimensionless diffusion equation:

$$\frac{\partial c}{\partial t} = \Delta c $$
with $\Delta$ the Laplace operator.

## Vertical diffusion:
We start by testing the vertical diffusion:

The initial condition is a step function:

~~~gnuplot
set multiplot layout 1,2
set xlabel 'x'
set ylabel 'c'
p[-0.5:0.5][-0.1:1.1] 'initialcondition.txt' u 1:3 t 'horizontal'
set xlabel 'y'
set ylabel 'c'
p[0.0:1.0][-0.1:1.1] 'initialcondition.txt' u 2:4 t 'vertical'
unset multiplot
~~~

Self-similarity:

~~~gnuplot
reset

set multiplot layout 2,1
set xlabel '\eta_x'
set ylabel 'c'
p "bulkH.txt" u ($2/sqrt($1)):3 t 'num',\
  0.5*erfc(x/2) t 'theo'
set xlabel '\eta_y'
set ylabel 'c'
p "bulkV.txt" u (($2-0.5)/sqrt($1)):3 t 'num',\
  0.5*erfc(x/2) t 'theo'
unset multiplot
~~~

# Code
We use the [multilayer hydrostatic solver](basilisk.fr/src/layered/hydro.h) in
1D. The `diffusion_neumann.h` file is equivalent to the
[diffusion.h](http://basilisk.fr/src/layered/hydro.h) file but with a different
boundary condition at the bottom.


## Libraries and Definitions
*/

//Librairies
#include "grid/multigrid1D.h"
#include "layered/hydro.h"

#include "./solvers/diffusion_neumann.h"


//Physical parameters
#define HINIT 1.
#define DIFFUSIVITY 1.
#define TMAX .1

//usefull variables and files
double altitude;
FILE * ic, * bulkV, * bulkH;

char name_init[30];
char name_bulkV[30];
char name_bulkH[30];

/**
## Main function
*/
int main(){
  L0 = 1.;
  X0 = -0.5;
  G=1.;
  nl = 64;
  N = 128;
  DT = 1e-5; // To be set manually, for stability of horizontal diffusion

  sprintf(name_init,"initialcondition.txt");
  ic = fopen(name_init,"w");
  fclose(ic);
  sprintf(name_bulkV,"bulkV.txt");
  bulkV = fopen(name_bulkV,"w");
  fclose(bulkV);
  sprintf(name_bulkH,"bulkH.txt");
  bulkH = fopen(name_bulkH,"w");
  fclose(bulkH);

  run();

}

/**
## Definition of the concentration field in the layered framework
*/
scalar c, c1;

event defaults(i=0){
  c1 = new scalar[nl];
  c  = new scalar[nl];
  tracers = list_append (tracers, c);
  tracers = list_append (tracers, c1);
}

event cleanup(t= end, last){
  delete({c,c1});
}

/**
## Diffusion step: 
Definition of the values of the fluxes at the top and the bottom of the domain
 */
double dcb=0., dct= 0.;

event diffusionOfC (i++){
  dt = dtnext(DT);
  horizontal_diffusion({c,c1},DIFFUSIVITY, dt);
  foreach(){
    vertical_diffusion_neumann(point, h, c , dt, DIFFUSIVITY,
                               dct, dcb);
    vertical_diffusion_neumann(point, h, c1, dt, DIFFUSIVITY,
                               dct, dcb);
  }
}

/**
## Initial condition: step function
*/

event init(i=0){
  ic = fopen(name_init,"a");

  foreach(){
    zb[]=0.;
    altitude =0.;
    foreach_layer(){
      h[] = HINIT/nl;
      altitude += h[]/2;
      c1[] = (altitude < HINIT/2.) ;
      c[] = (x < 0);
      fprintf(ic,"%g %g %g %g\n",x, altitude, c[], c1[]);
      altitude += h[]/2;

    }
  }
  fclose(ic);
}
  
/**
## Outputs
*/
event output_data(t += .001, t<= TMAX){
  bulkV = fopen(name_bulkV,"a");
  bulkH = fopen(name_bulkH,"a");

  foreach(){
    altitude =0.;
    foreach_layer(){
      altitude += h[]/2;
      if(x>-0.1 && x< 0.1 )
        fprintf(bulkV,"%g %g %g\n",t, altitude, c1[]);
      if(altitude >0.4 && altitude < 0.6)
        fprintf(bulkH,"%g %g %g\n",t, x, c[]);
      altitude += h[]/2;
    }
  }
  fclose(bulkV);
  fclose(bulkH);
}

event output_progress(i += 50){
  fprintf(stderr, "i=%d, t= %g\n",i,t);
}
    