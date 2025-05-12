/**
# Advection of a Gaussian Pulse on an Equidistant Grid
On this page, a check how spatial discretization introduces numerical errors for a simple test case concerning the advection solver is presented. The case will consist of a one-dimensional Gaussian pulse that advects from $x=-0.5$ to $x=0.5$ over time. The case is setup very similar to the diffusion test presented [here](diffgausspulse.c) 
*/
#include "grid/bitree.h"
#include "run.h"
#include "tracer.h"

double t0=1./3.;
double xo=-0.5;
face vector uf[];
scalar c[];
scalar * tracers = {c};
int m;
char name1[100];
FILE * fp1;

double f(double t,double to,double xi,double Delta){
  return (pow(M_PI*0.01*to,0.5)*(erf(((xi-xo-t)+(Delta/2))/(2*pow(0.01*(to),0.5)))-erf(((xi-xo-t)-(Delta/2))/(2*pow(0.01*(to),0.5)))))/Delta;
}
/**
Seven simulations are run with increasing grid resolution.
*/
int main(){
  sprintf(name1,"errvsD.dat");
  fp1 = fopen(name1,"w");
  L0=5;
  X0=-L0/2;
  for (m=0;m<7;m++)
    run(); 
}

event init(t=0){
  DT=1e-4; 
  init_grid(1<<(6+m));
  foreach()
    c[]=f(t,t0,x,Delta);
  foreach_face()
    uf.x[]=1.;
}
/**
Advection is solved for by default and does not require an user-defined event. We do set the timestep each iteration in a new event.
*/
event timestepper(i++;t<=1.) {
  dt=dtnext(DT);
}

event end(t=1){
  char name[100];
  sprintf(name,"advGauss%d.dat",m);
  FILE * fp = fopen(name,"w");
  foreach()
    fprintf(fp,"%g\t%g\t%g\t%g\n",x,c[],f(t,t0,x,Delta),fabs(c[]-f(t,t0,x,Delta)));
  fflush(fp);
  double toterr=0;
  foreach()
    toterr+=fabs(c[]-f(t,t0,x,Delta))*Delta;
  fprintf(fp1,"%d\t%g\t%g\n",m,(double) L0/pow(2.,6.+(double)m),toterr);
  fflush(fp1);
}
/**
## Results
We check the obtained solutions at $t=t_{\mathrm{end}}$. 

 ~~~gnuplot
  set xr [-1:1]
  set yr [-0.1:1.4]
  set xlabel 'x'
  set ylabel 'c'
  set size square
  set samples 1000
    plot exp(-(x-0.5)**2/(4*0.01*1/3)) with line lw 5 lc rgb 'magenta' title "Analytical solution" ,\
         'advGauss1.dat' using 1:2 with line lw 3 lc rgb 'red'  title "N=128" ,\
         'advGauss2.dat' using 1:2 with line lw 3 lc rgb 'blue' title  "N=256" ,\
         'advGauss3.dat' using 1:2 with line lw 3 lc rgb 'green' title  "N=512" ,\
         exp(-(x+0.5)**2/(4*0.01*1/3)) lc rgb 'black' title 'Initialized solution'
  ~~~

The solutions do not look very accurate for the coarser-grid runs. Lets study the erros:

 ~~~gnuplot
  set xr [-0.1:1.1]
  set yr [0:0.25]
  set xlabel 'x'
  set ylabel 'Error'
  set size square
    plot 'advGauss1.dat' using 1:4 with line lw 3 title "N=128 run" ,\
         'advGauss2.dat' using 1:4 with line lw 3 title "N=256 run" ,\
         'advGauss3.dat' using 1:4 with line lw 3 title "N=512 run" 
   ~~~
  
We see that increasing the grid resolution caused a drastic decrease of the error, eventough the overal error-distribution structure seems similar. Nice convergence properties, lets study these in a bit more detail:

~~~gnuplot
  set xr [0.001:0.1]
  set yr [0.00002:0.5]
  set logscale y
  set logscale x
  set xlabel 'D'
  set ylabel 'Total Error'
  set key left top box 3
  set size square
    plot    (60*x**2) lw 3 lc rgb 'purple' title 'second order}',\
             'errvsD.dat' using 2:3 pt 4 title 'Error'
~~~

We observe a well behaved second-order accuracy in space. Well done src/tracer.h! Notice that the coarser-grid results lays outside the so-called 'convergence region' of the solver.

## The next step
The next step is to use adaptive grid refinement for this case. We will check if the leassons learned with the diffusion solver are transferable. 
*/