/**
#Convergence test forthe ekman spiral using an adaptive grid
Very similar to the approach we took for the [fixed-grid](ekman.c) convergence testing, we check the spatial order of convergence of the ‘dynamical-core’ of the single collumn model, i.e. the diffusion solver. We use a one-dimensional bi-tree adaptive grid and the generic timeloop header file together with the aforementioned solver.
*/
#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"
/**
For reasons to do with brevity, this page only accentuates the differences in the set-up compared to the fixed-grid-testing approach. If you have any questions regarding the set-up and you have not read the corresponding fixed-grid-test page, than please look [there](ekmanfg.c) first.  
*/
double ekmanexv (double x){
  return sin(x)*exp(-x);
}

double ekmanexu (double x){
  return 1-(cos(x)*exp(-x));
}

double Ekmanvdx (double x){
  return -0.5*exp(-x)*(cos(x)+sin(x));
}

double Ekmanudx (double x){
  return 0.5*exp(-x)*(2.*exp(x)*x-sin(x)+cos(x));
}

double solv(double x, double delt){
  double x1 = x-delt/2.;
  double x2 = x+delt/2.; 
  return (Ekmanvdx(x2)-Ekmanvdx(x1))/delt;
}

double solu(double x, double delt){
  double x1=x-delt/2.;
  double x2=x+delt/2.;
  return (Ekmanudx(x2)-Ekmanudx(x1))/delt;
}


int maxlevel;
scalar u[],v[];
const face vector mum[]={0.5};
int j;
FILE * fpo;
/**
The used adaptation algorithm requires an *error threshold* or *refinenent criterion* for the velocity components $\zeta$, that is declared below.
*/
double zeta;

int main()	
{
  X0=0.;
  L0=100.;
  fpo = fopen ("overviewagekman","w");
  /**
  The calculations on this case are repeated 20 times over so that the maximum grid size and refinement criterion can be reduced iteratively. 
  */
  for (j=0;j<20;j++){
    run();
  }
}	

event init(t=0){
  u[left]=dirichlet(0.);
  v[left]=dirichlet(0.);
  u[right]=dirichlet(ekmanexu(L0));
  v[right]=dirichlet(ekmanexv(L0));
  /**
  The grid is initialized with a relatively coarse resolution, employing 32 points. The maximum level is doubled each iteration and the refinement criterion is halved each iteration, starting from 0.1. 
  */
  init_grid(1<<(4));
  maxlevel = 5+j;
  zeta = 0.1/pow(2,(double)j);
  dt=0.01;
  DT=0.01;
  foreach(){
    u[]=solu(x,Delta);
    v[]=solv(x,Delta);
  }
  boundary(all);
  /**
  Before the run is started, the solution is initialized on a grid that is consistent with the adaptations settings. 
  */
  while(adapt_wavelet({u,v},(double []){zeta,zeta},maxlevel).nf){
    foreach(){
      u[]=solu(x,Delta);
      v[]=solv(x,Delta);
    }
    boundary(all);
  }
}

event Diffusion(i++){
  scalar rx[],ry[];
  foreach(){
    rx[]=v[];
    ry[]=(1.-u[]);
  }
  boundary({ry,rx});
  dt=dtnext(DT);
  diffusion(u,dt,mum,rx);
  diffusion(v,dt,mum,ry);
}
/**
The fidelity in the representation of the solution by the grid is evaluated each timestep. The grid is adapted, if necessarry.  
*/
event adapt(i++){
  adapt_wavelet({u,v},(double []){zeta,zeta},maxlevel);
}
/**
The details of the solution are evaluated at the end of each `refinement iteration' and written to a file named; *ekmanadaptive.dat*. 
*/
event eval(t=(10)){
  static FILE * fp = fopen("ekmanadaptive.dat","w");
  double eta = 0;
  int n = 0;
  foreach(){
    eta+=fabs(u[]-solu(x,Delta))*Delta;
    eta+=fabs(v[]-solv(x,Delta))*Delta;
    n++;
  }
  fprintf(fp,"%d\t%g\t%g\n",n,eta,perf.t);
  fflush(fp);
  scalar wu[], wv[];
  wavelet(u,wu);
  wavelet(v,wv);
  foreach()
    fprintf(fpo,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n",x,Delta,fabs(u[]-solu(x,Delta)),fabs(v[]-solv(x,Delta)),wu[],wv[],zeta);
  fflush(fpo);
}
/**
## Results
We can check if everything has gone according to the intended set-up by plotting the number of grid cells and the corresponding error for each iteration.

~~~gnuplot
set logscale x
set logscale y
set xlabel 'Number of grid cells'
set ylabel 'Total error'
plot 'ekmanadaptive.dat' u 1:2 w l lw 3 t 'Adaptive' ,\
      x**-2 lw 3 t '{/Symbol \265} N^{-2}'
~~~

We can conclude that the error (agian) scales inversely propotional to the square of the number of grid cells. We also look at the error vs esitmated error in `v`:

~~~gnuplot this is (in part) the same data as fig. 2a in the GMD paper
set xlabel 'Xsi for v'
set ylabel 'Error in v'
set key off
plot 'overviewagekman' u 4:6 pt 5
~~~

We observe a correlation along the 1:1 line. Compared to the plot in the paper, the are additonal scatter points with very low error, I doubt this is due to the newer version of Basilisk. I will look into this...
*/