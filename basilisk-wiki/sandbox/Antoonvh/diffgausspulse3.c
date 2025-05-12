/**
# Diffusion of a Gaussian Pulse Using an Adaptive Grid
As a follow-up of [this](diffgausspulse2.c) page, we study the effect of using adaptive-grid refinement on the accuracy of the solution. On the aforementioned page we found that a resolution boundary may introduce numerical errors due to the prolongation at resolution boundaries. We will see that an adaptive-grid approach is not free its own additional issues. Rather than being directly controlled by the initialized grid resolution, the accuracy of the solution is now a function of the so-called refinement criterion ($\zeta$). This criterion will in turn dictate the (spatiotemporal varying) grid resolution.
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
We run four different simulations. We always use $\zeta=0.05$, unless stated otherwise. Supplemented with a maximum resolution corresponding to $N=128$.

1. Using a fixed grid with $N=128$ 
2. Using the Default settings 
3. Using $\zeta$=0.025
4. Using $refine\_linear$ for prolongation
 
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
  if (m>0)
    zeta=0.05;
  if (m==2)
    zeta=0.025;
  if (m>2)
    c.prolongation=refine_linear;
  foreach()
    c[]=f(t,t0,x,D,Delta);
  if (m>0){
    boundary({c});
    while(adapt_wavelet({c},(double[]){zeta},7,5).nc){
      foreach()
        c[]=f(t,t0,x,D,Delta);
      boundary({c});
    }
  }
}

event diffn(i++;t<=1.) {
  dt=dtnext(DT);
  diffusion(c,dt,DD); 
}
/**
An adapt event is added for adaptivity purposes.
*/
event adapt(i++){
  if (m>0){
    adapt_wavelet({c},(double[]){zeta},7,5); 
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
  plot 'Gauss0.dat' using 1:2 ps 3 lc rgb 'red'  title "Fixed-grid run" ,\
       'Gauss1.dat' using 1:2 ps 2 lc rgb 'blue' title  "{/Symbol z}=0.05" ,\
       'Gauss2.dat' using 1:2 ps 1 lc rgb 'green' title "{/Symbol z}=0.025" 
 ~~~
  
Looks as can be expected. Lets look at the solution from the first two adaptive-grid runs.

 ~~~gnuplot
  unset logscale y
  set xr [-1 :1]
  set yr [0 : 1.1] 
  
  set xlabel 'x'
  set ylabel 'c'
  set size square
  set samples 10001
    plot  0.5*exp(-((x)**2)/(0.04*4/3)) lw 5 lc rgb 'magenta' title "Exact solution" ,\
         'Gauss1.dat' using 1:3 with line lw 3 lc rgb 'blue' title "{/Symbol z}=0.05}" ,\
         'Gauss2.dat' using 1:3 with line lw 2 lc rgb 'green' title "{/Symbol z}=0.025}"  ,\
         exp(-((x)**2)/(0.04*1/3)) lc rgb 'black' with line title "Initialized pulse"
~~~
  
This looks O.K.-ish. Lets study the errors for all runs in some more detail. 

 ~~~gnuplot
  set xr [-1:1]
  set yr [0:0.01]
  set xlabel 'x'
  set ylabel 'Error'
  set size square
    plot 'Gauss0.dat' using 1:5 with line lw 3 lc rgb 'red' title "Fixed-grid run" ,\
         'Gauss1.dat' using 1:5 with line lw 3 lc rgb 'blue' title "{/Symbol z}= 0.05 "  ,\
         'Gauss2.dat' using 1:5 with line lw 2 lc rgb 'green' title "{/Symbol z}= 0.025" 
~~~

This plot reveals some interesting results. First, the fixed-grid run is most accurate compared to the adaptive-grid counterparts. Second, when the $\zeta$ is set to a smaller value, the obtained error do not really appear to decrease.

## Analysis: Why does using a smaller $\zeta$-value not yield better results?
From the first figure it can be seen that setting a lower refinement cirterion does lead to more refinement. We have learned from [this](diffgausspulse.c) page that higher resolution grids can be associated with a more accurate solution. However, [this](diffgausspulse2.c) page presented a cautionary note on the introduction of resolution boundaries. For this case, with the chosen $\zeta$-values, the increased accuracy of using more cells is not setoff by the error introduced with the additional resolution boundaries. With the aim to check this, a simulation was run with using $refine\_linear$ prolongation attribute. The results are shown below:

 ~~~gnuplot
  set xr [-1:1]
  set yr [0:0.01]
  set xlabel 'x'
  set ylabel 'Error'
  set size square
   plot   'Gauss0.dat' using 1:5 with line lw 3 lc rgb 'red' title "Fixed-grid run" ,\
         'Gauss1.dat' using 1:5 with line lw 3 lc rgb 'blue' title "{/Symbol z}= 0.05 "  ,\
         'Gauss3.dat' using 1:5 with line lw 2 lc rgb 'orange' title "linear prolongation,{/Symbol z}= 0.05" 
~~~

It seems to be doing better, but not too convincingly. Lets check the used grid for both simulations at $t=t_{end}$.

~~~gnuplot
set xr [-1:1]
set yr [0.03:0.4]
  set xlabel 'x'
  set ylabel '{/Symbol D}'
  set logscale y
  set size square
  plot 'Gauss1.dat' using 1:2 ps 3 lc rgb 'blue'  title "{/Symbol z}=0.05" ,\
       'Gauss3.dat' using 1:2 ps 2 lc rgb 'orange' title  "linear prolongation,{/Symbol z}=0.05" 
~~~

Obviously, the linear-prolongation run is using much fewer cells. The is the result of the fact that linear prolongation is more accurate! The corresponding wavelet estimated error will in general be smaller compared to the (default) bilinear prolongation. This in turn will lead to more coarsening and less refinement for the more accurate method. Therefore, the comparison of the errors made in this section is not really fair. 

## The next step
To circumvent this issue a work arround can be introcuded such that $refine\_bilinear$ is used when making the wavelet-based grid assessment whilst using $refine\_linear$ for the treatment of resolution boundaries. The corresponding setup and results can be found [here](diffgausspulse3b.c). It was found that the presented results here, represent a rather special case.   
*/