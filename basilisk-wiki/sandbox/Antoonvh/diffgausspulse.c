/**
# Diffusion of a Gaussian Pulse on an Equidistant Grid
On this page, a check how spatial discretization introduces numerical errors for a simple test case concerning the diffusion solver is presented. The case will consist of a one-dimensional Gaussian pulse that diffuses over time.
*/
#include "grid/bitree.h" //The tree-based grid will come in handy later
#include "run.h"
#include "diffusion.h"
/**
The case is setup such that the analytical solution will be,
$$c(x,t)=\sqrt{\frac{t_0}{t+t_0}}e^{-\frac{x^2}{4D(t+t_0)}},$$
with $t_0$ and $D$ some time and diffusivity related constants, respectively. Since the discretized solution ($c_i$) will represent volume averaged values, we need a function that translates the analytical solution to a locally averaged one. 
  
$$c_i=\frac{1}{\Delta}\sqrt{\frac{t_0}{t+t_0}}\int_{x_i-\Delta/2}^{x_i+\Delta/2}e^{-\frac{x^2}{4D(t+t_0)}}\mathrm{d}x=$$
$$\frac{1}{\Delta}\sqrt{t_0D\pi} \left[\mathrm{erf}\left(\frac{x}{2\sqrt{D(t+t_0)}}\right)\right]_{x_i-\Delta/2}^{x_i+\Delta/2},$$
with $x_i$ the $x$-coordinate of the cell centre, $\Delta$ the corresponding grid box size and $\mathrm{erf()}$ the so-called error funtion.  
*/
double f(double t,double to,double xi,double D,double Delta){
  return (pow((to)/(to+t),0.5)*pow(M_PI*D*(t+to),0.5)*(erf((xi+(Delta/2))/(2*pow(D*(t+to),0.5)))-erf((xi-(Delta/2))/(2*pow(D*(t+to),0.5)))))/Delta;
}
/**
This function will be used to initialize the numerical solution and evaluate the error at the beginning and end of the simulation, respectively. Now we arrive at the simulation itself. First, some handy constants are defined and fields are allocated.  
*/
double t0=1./3.;
double D=0.01;
const face vector DD[]= {D};
scalar c[];
int m;
char name1[100];
FILE * fp1;
/**
Second, the domain is setup and the simulation is initialized with different resolutions, varying from $N=64$ to $N=4096$.
*/
int main(){
  sprintf(name1,"toterrvsD.dat");
  fp1=fopen(name1,"w");
  L0=5;
  X0=-L0/2;
  for (m=0;m<7;m++)
    run(); 
}

event init(t=0){
  DT=1e-4; //Erros due to time integration could be limited. 
  init_grid(1<<(6+m));
  foreach()
    c[]=f(t,t0,x,D,Delta);
}
/**
Third, Diffusion is solved for, such that by construct of the chosen constants, the height of the Gaussian pulse in half its initial value at the end of the simulation. 
*/
event diffn(i++;t<=1.) {
  dt=dtnext(DT);
  diffusion(c,dt,DD); 
}
/**
Finally, A file is written with the solution and the corresponding error, such that we can analyze it once the simulation is finished. 
*/
event end(t=1){
  char name[100];
  sprintf(name,"Gauss%d.dat",m);
  FILE * fp = fopen(name,"w");
  foreach()
    fprintf(fp,"%g\t%g\t%g\t%g\n",x,c[],f(t,t0,x,D,Delta),fabs(c[]-f(t,t0,x,D,Delta)));
  fflush(fp);
  double toterr=0;
  foreach()
    toterr+=fabs(c[]-f(t,t0,x,D,Delta))*Delta;
  fprintf(fp1,"%d\t%g\t%g\n",m,L0/((double)pow(2.,(double)m+6.)),toterr);
  fflush(fp1);
}
/**
## Results
We check the obtained solutions at $t=t_{\mathrm{end}}$. 

~~~gnuplot
  set xr [-1:1]
  set yr [0:1.1]
  set xlabel 'x'
  set ylabel 'c'
  set size square
  set samples 1000
    plot 'Gauss0.dat' using 1:2 with line lw 6  title "Coarse-grid run solution" ,\
         'Gauss1.dat' using 1:2 with line lw 4 title  "Fine-grid run solution" ,\
         0.5*exp(-x**2/(4*0.01*4/3)) with line lw 2 title "Analytical solution" ,\
         exp(-x**2/(4*0.01*1/3)) lc rgb 'black' title 'Initialized solution'
~~~
  
Looks O.K. Lets zoom in on the error,

~~~gnuplot
  set xr [-1:1]
  set yr [0:0.01]
  set xlabel 'x'
  set ylabel 'Error'
  set size square
    plot 'Gauss0.dat' using 1:4 with line lw 3 title "Coarse-grid run error" ,\
         'Gauss1.dat' using 1:4 with line lw 3 title "Fine-grid run error" 
~~~
  
We see that increasing the grid resolution caused a drastic decrease of the error, eventough the overal error-distribution structure seems similar. Lets check these convergence properties.

~~~gnuplot
  set xr [0.001:0.1]
  set yr [0.000002:0.01]
  set logscale y
  set logscale x
  set xlabel '{/Symbol D}'
  set ylabel 'Total Error'
  set key left top box 1
  set size square
   plot    (0.5*x**2) lw 3 lc rgb 'purple' title '{/Symbol \265}{/Symbol D}^{2}',\
            'toterrvsD.dat' using 2:3 pt 4 title 'Total error' 
~~~
         
We observe second-order spatial accuracy. Very nicely behaved, well done $diffusion.h$! The run with $N=4096$ does not seem to follow the purple scaling line. There might be two reasons (in order of likelyhood):

1. The errors due to the time integration scheme start to overwelm the errors in the solution.
2. The precision of the machine is limiting the accurate evaluation on the solution. 

## The next step
It appears that the error is the solution is only appreciable in a subsection of the domain. Furthermore, the centre of the Gaussian curve seems to benefit most from the increased resolution. Hence, it seems sensible to locally refine the grid there. Fortunately, the tree-based grid facilitates a convinient implementation of anisotropic grids. The definition and the results of this next step can be found [here](diffgausspulse2.c).
*/