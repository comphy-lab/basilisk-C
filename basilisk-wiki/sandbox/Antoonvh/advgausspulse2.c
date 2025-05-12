/**
# Advection of a Gaussian Pulse on an Adaptive Grid
Similar to a study for the diffusion solver [here](diffgausspulse3b.c), this page presents how adaptivity affects the error in the solution for an advection problem (introduced [here](advgausspulse.c). The focus will be on the role of the refinement criterion $\zeta$. 
*/
#include "grid/bitree.h"
#include "run.h"
#include "tracer.h"

double t0=1./3.;
double xo=-0.5;
int maxlevel = 12;
double zeta,referr;
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
Fourteen simulations are run with different values for the refinement criterion $\zeta$, supplemented with a benchmark run using a fixed-equidistant grid with $N=4096$. 

*/
int main(){
  L0=5;
  X0=-L0/2;
  sprintf(name1,"zetavstoterr.dat");
  fp1=fopen(name1,"w");
  for (m=0;m<15;m++)
    run(); 
}

event init(t=0){
  DT=1e-4; 
  init_grid(1<<maxlevel);
  zeta=0.1/((double) pow(2.,(double)m));
  if (m==0)
    zeta=0.;
  foreach()
    c[]=f(t,t0,x,Delta);
  while (adapt_wavelet({c},(double[]){zeta},maxlevel,2,{c}).nc){
    foreach()
      c[]=f(t,t0,x,Delta);
    boundary({c});
  }
  foreach_face()
    uf.x[]=1.;
}

event timestepper(i++;t<=1.) {
  dt=dtnext(DT);
}

event adapt(i++){
  adapt_wavelet({c},(double[]){zeta},maxlevel,2,{c});
  foreach_face()
    uf.x[]=1.; 
}
/**
The total error is also evaluated and the difference with the benchmark result is written to a file.
*/
event end(t=1){
  char name[100];
  sprintf(name,"aadvGauss%d.dat",m);
  FILE * fp = fopen(name,"w");
  foreach()
    fprintf(fp,"%g\t%g\t%g\t%g\n",x,c[],f(t,t0,x,Delta),fabs(c[]-f(t,t0,x,Delta)));
  fflush(fp);
  double toterr=0;
  foreach()
    toterr+=fabs(c[]-f(t,t0,x,Delta))*Delta;
  if (m==0)
    referr=toterr;
  else{
    fprintf(fp1,"%d\t%g\t%g\n",m,zeta,toterr-referr);
    fflush(fp1);
  }
}
/**
## Results
We check the obtained errors at $t=t_{\mathrm{end}}$. 

 ~~~gnuplot
  set xr [-0.1:1.1]
  set yr [0:0.11]
  set xlabel 'x'
  set ylabel 'Error'
  set size square
    plot 'aadvGauss2.dat' using 1:4 with line lw 3 title "{/Symbol z} = 0.025" ,\
         'aadvGauss3.dat' using 1:4 with line lw 3 title "{/Symbol z} = 0.0125" ,\
         'aadvGauss4.dat' using 1:4 with line lw 3 title "{/Symbol z} = 0.00625" ,\
         'aadvGauss0.dat' using 1:4 with line lw 2 lc rgb 'black' title "Fine-grid run error"
         ~~~
   
We see that decreasing the refinement criterion leads to a decrease of the error. More specifically, the results seem to convergence towards the fixed-equidistant grid run results. Lets look at this convergence. 
~~~gnuplot
  
  set xr [0.000005:0.06]
  set yr [0.000001:0.3]
  set logscale y
  set logscale x
  set xlabel '{/Symbol z}'
  set ylabel 'Total Error - Benchmark error'
  set size square
    plot 'zetavstoterr.dat' using 2:3 title 'Total error - Benchmark error' ,\
    (2*x**1.2) lw 3 title '{/Symbol \265}{/Symbol z}^{1.2}'
         ~~~
         
We observe that the error converges towards the error that was obtained with the fixed-equidistant grid run In a nicely behaved manner. Well done $adapt\_wavelet$ in conjunction with src/tracer.h! 

## The next step

*/