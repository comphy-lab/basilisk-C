/**
# Diffusion of a Gaussian Pulse Using an Adaptive Grid B
As a side note on [this](diffgausspulse3.c) page, we study the convergence of the error in the solution with smaller refinement criterion ($\zeta$) values. Notice that we will take a equidistant-fixed-grid run as reference to compare our errors against.
*/
#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

double t0=1./3.;
double D=0.01;
double zeta,referr;
const face vector DD[]= {D};
scalar c[];
int m,mm;
int maxlevel = 11;
FILE * fp1;
char name1[100];

double f(double t,double to,double xi,double D,double Delta){
  return (pow((to)/(to+t),0.5)*pow(M_PI*D*(t+to),0.5)*(erf(((xi)+(Delta/2))/(2*pow(D*(t+to),0.5)))-erf(((xi)-(Delta/2))/(2*pow(D*(t+to),0.5)))))/Delta;
}
/**
We run nine different simulations with varying $\zeta$, supplemented with a benchmark run using the maximum resolution corresponding to $N=2048$.
*/
int main(){
  L0=5;
  X0=-L0/2;
  for (mm=0;mm<2;mm++){
    sprintf(name1,"zetavstoterr%d.dat",mm);
    if (mm==1)
      fclose(fp1);
    fp1=fopen(name1,"w");
    for (m=0;m<10;m++)
      run();
  }
}

event init(t=0){
  DT=1e-4;
  init_grid(1<<maxlevel);
  zeta=0.1/((double) pow(2.,(double)m));
  if (m==0)
    zeta=0.0;
  foreach()
    c[]=f(t,t0,x,D,Delta);
  if (m>0){
    boundary({c});
    while(adapt_wavelet({c},(double[]){zeta},maxlevel,2).nc){
      foreach()
        c[]=f(t,t0,x,D,Delta);
      boundary({c});
    }
    if (mm){
      c.prolongation=refine_linear;
      c.refine=refine_linear;
      boundary({c});
    }
  }
}

event diffn(i++;t<=1.) {
  dt=dtnext(DT);
  diffusion(c,dt,DD); 
}

event adapt(i++){
  if (m>0){
    if (mm){
      c.prolongation=refine_bilinear; //Back to default
      boundary({c});
    }
    adapt_wavelet({c},(double[]){zeta},maxlevel,2); 
    if (mm){
      c.prolongation=refine_linear; //Undo the default setting again
      boundary({c});
    }
  }
}

event end(t=1){
  char name[100];
  sprintf(name,"Gauss%d.dat",m);
  FILE * fp = fopen(name,"w");
  foreach()
    fprintf(fp,"%g\t%g\t%g\t%g\t%g\n",x,Delta,c[],f(t,t0,x,D,Delta),fabs(c[]-f(t,t0,x,D,Delta)));
  fflush(fp);
  double toterr=0;
  foreach()
    toterr+=fabs(c[]-f(t,t0,x,D,Delta))*Delta;
  if (m==0)
    referr=toterr;
  else{
    fprintf(fp1,"%d\t%g\t%g\n",m,zeta,toterr-referr);
    fflush(fp1);
  }
}
/**
## Results
Lets study how the total error compared to the benchmark result scales with $\zeta$.  

~~~gnuplot
  
  set xr [0.0001:0.1]
  set yr [0.000002:0.01]
  set logscale y
  set logscale x
  set xlabel '{/Symbol z}'
  set ylabel 'Total Error - Benchmark error'
  set key left top box 3
  set size square
    plot    (0.05*x**0.7) lw 3 lc rgb 'purple' title '{/Symbol \265}{/Symbol z}^{0.7}',\
            (0.05*x**1.0) lw 3 lc rgb 'green' title '{/Symbol \265}{/Symbol z}^{1.0}',\
            'zetavstoterr0.dat' using 2:3 pt 4 title 'Default settings' ,\
            'zetavstoterr1.dat' using 2:3 pt 4 title 'With linear attributes' 
   
         ~~~
         
The absolute error obtained with the non-default settings is always lower than those obtained with the default ones. Furthermore, the scaling is along a flatter line for the latter. Notice that the scaling is not well behaved since the results do not line-up properly allong their approximate scaling order.  
*/