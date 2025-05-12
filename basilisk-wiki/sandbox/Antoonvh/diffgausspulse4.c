/**
# Diffusion of a Gaussian Pulse Using an Adaptive Grid, revisited.
As a follow-up of [this](diffgausspulse3.c) page, we study the effect of using linear prolongation for the treatment of resolution boundaries when using an adaptive grid. We have learned:

* Higher resolutions yield more accurate results.
* Linear prolongation yields more accurate results when treating resolution boundaries. $O(\Delta^{3})$-accuracy in 1D.
* One should be cautionary when setting up a case with adaptive refinement on the excessive introduction of resolution boundaries.

*/
#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

double t0=1./3.;
double D=0.01;
double zeta;
const face vector DD[]= {D};
scalar c[];
int m;

double f(double t,double to,double xi,double D,double Delta){
  return (pow((to)/(to+t),0.5)*pow(M_PI*D*(t+to),0.5)*(erf(((xi)+(Delta/2))/(2*pow(D*(t+to),0.5)))-erf(((xi)-(Delta/2))/(2*pow(D*(t+to),0.5)))))/Delta;
}
/**
We run three different simulations, using $\zeta=0.05$, unless stated otherwise.

1. Using the Default settings
2. Using linear prolongation for resolution boundaries, but bilinear interpolation when making the wavelet-based grid assessment 
3. Similar to 2, but now with $\zeta$=0.025
4. Similar to 3, but with $refine\_linear$ for the refinement attribute
 
*/
int main(){
  L0=5;
  X0=-L0/2;
  for (m=0;m<4;m++)
    run(); 
}

event init(t=0){
  
  DT=1e-4;
  init_grid(1<<7);
  zeta=0.05;
  if (m>=2)
    zeta=0.025;
  foreach()
    c[]=f(t,t0,x,D,Delta);
  while(adapt_wavelet({c},(double[]){zeta},7,5).nc){
    foreach()
      c[]=f(t,t0,x,D,Delta);
    boundary({c});
  }
  if(m>0){
    c.prolongation=refine_linear;
    boundary({c});
  }
  if(m==3)
    c.refine=refine_linear;
}

event diffn(i++;t<=1.) {
  dt=dtnext(DT);
  diffusion(c,dt,DD); 
}

event adapt(i++){
  if(m>0){
    c.prolongation=refine_bilinear;
    boundary({c});
  }
    adapt_wavelet({c},(double[]){zeta},7,5); 
  if(m>0){
    c.prolongation=refine_linear;
    boundary({c});
  }
}

event end(t=1){
  char name[100];
  sprintf(name,"Gauss%d.dat",m);
  FILE * fp = fopen(name,"w");
  foreach()
    fprintf(fp,"%g\t%g\t%g\t%g\t%g\n",x,Delta,c[],f(t,t0,x,D,Delta),fabs(c[]-f(t,t0,x,D,Delta)));
  fflush(fp);
}
/**
## Results
First, the grid resolution at $t=t_{end}$ for some different runs are plotted. 

~~~gnuplot
set xr [-1:1]
set yr [0.03:0.4]
  set xlabel 'x'
  set ylabel '{/Symbol D}'
  set logscale y
  set size square
  plot 'Gauss0.dat' using 1:2 ps 3 lc rgb 'red'  title "Default" ,\
       'Gauss1.dat' using 1:2 ps 2 lc rgb 'blue' title  "Linear treatment of res. bound." ,\
       'Gauss2.dat' using 1:2 ps 1 lc rgb 'green' title "Idem. {/Symbol z}=0.025" 
 ~~~
  
Now the grids corresonding to the runs woth $\zeta=0.05$ are similar. Lets study the errors for all runs in some more detail. 

 ~~~gnuplot
  unset logscale y
  set xr [-1:1]
  set yr [0:0.01]
  set xlabel 'x'
  set ylabel 'Error'
  set size square
    plot 'Gauss0.dat' using 1:5 with line lw 3 lc rgb 'red' title "Default" ,\
         'Gauss1.dat' using 1:5 with line lw 3 lc rgb 'blue' title "Linear treatment of res. bound."  ,\
         'Gauss2.dat' using 1:5 with line lw 3 lc rgb 'green' title "Idem. {/Symbol z}= 0.025" 
~~~

We now see that the linear prolongation attribute for resolution boundaries works quite well, also on adaptive grids. However, decreasing refinement cirterion to $\zeta=0.025$, did not yield much better results. Also there might be some cancelling out errors with opposite sign going on. Furhtermore, as the pulse widened over the course of the simulation with $\zeta=0.025$, the grid was also locally refined. This means that the $c.refine$ attribute was used at some point. We compare to results for the default settings and with $refine\_linear$ for refinement.  

 ~~~gnuplot
  set xr [-1:1]
  set yr [0:0.01]
  set xlabel 'x'
  set ylabel 'Error'
  set size square
  plot  'Gauss2.dat' using 1:5 with line lw 3 lc rgb 'green' title "Linear prolongation{/Symbol z}= 0.025" ,\
        'Gauss3.dat' using 1:5 with line lw 3 lc rgb 'orange' title "Linear refinement and prolongation {/Symbol z}= 0.025" 
~~~

This appears to be an improvement. Now the magnitude of the error corresponds somewhat with the fixed-grid run that we have performed [earlier](diffgausspulse.c).

The conslusion of this analysis is: More accurate numerical formulations result in more accurate solutions.   

## The next step
The next step is to study the advection of a Gaussian pulse. We will redo the analysis similar to as we have done for the diffusion solver, starting [here](404).  
*/
