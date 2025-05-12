/**
# Convergence test for the Ekman spiral using an equidistant grid. 
We check the spatial order of convergence of the 'dynamical-core' of the single collumn model, i.e. the diffusion solver. We use a one-dimensional Multigrid and the generic timeloop header file together with the aforementioned solver
*/
#include "grid/multigrid1D.h"
#include "run.h"
#include "diffusion.h"
/**
For initialization, setting boundary conditions and evaluating the accuracy of the end results we define functions that return the analytical solution. We use a normalized length scale $\gamma$ and geowind, so that, 

$$v=\mathrm{sin}(x)e^{-z},$$
$$u=1-\mathrm{cos}(z)e^{-z}.$$
*/

double ekmanexv (double z){
  return sin(z)*exp(-z);
}

double ekmanexu (double z){
  return 1-(cos(z)*exp(-z));
}
/**
Since the cell values in a finite-volume formulation represent length averages, we need to define the integral form, starting with the primitives:
*/
double Ekmanvdx (double x){
  return -0.5*exp(-x)*(cos(x)+sin(x));
}

double Ekmanudx (double x){
  return 0.5*exp(-x)*(2.*exp(x)*x-sin(x)+cos(x));
}
/**
And the integrals:
*/

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
/**
The scalar fields for the two wind components are declared here, together with the viscosity ($\nu$). Since we use a normalized $\gamma$ and Coriolis parameter we set $\nu=0.5$.  
*/
scalar u[],v[];
const face vector nu[]={0.5};
int j;
FILE * fpo;

/**
The simulation is run 10 times with increasing *j* values. We also set the boundary conditions consistent with a no-slip bottom (dubbed *left*) and top boundary conditions (dubbed *right*), consistent with the analytical solution.
*/
int main()	
{
  X0=0.;
  L0=100.;
  fpo = fopen("ekmanfgoverview", "w");
  for (j=0;j<10;j++){
    run();
  }
}	

u[left]=dirichlet(0.);
v[left]=dirichlet(0.);
u[right]=dirichlet(ekmanexu(L0));
v[right]=dirichlet(ekmanexv(L0));
/**
Where *j* controls the level of refinement, starting from a N=64 grid, iteratively doubling the number of grid cells to N= 32768. The timestep is set and the solution is initialized according to the Ekman solution.  
*/
event init(t=0){
  init_grid(1<<(j+6));
  dt=0.01;
  DT=0.01;
  foreach(){
    u[]=solu(x,Delta);
    v[]=solv(x,Delta);
  }
  boundary(all);
}
/**
Time integration is performed in this event. The Coriolis term is threated, in the formward direction, as a tendency term (*rx,ry*) for the Poisson-Helmholtz equations describing the evolution of *u* and *v*, respectively.  
*/
event Diffusion(i++){
  scalar rx[],ry[];
  foreach(){
    rx[]=v[];
    ry[]=(1.-u[]);
  }
  boundary({ry,rx});
  dt=dtnext(DT);
  diffusion(u,dt,nu,rx);
  diffusion(v,dt,nu,ry);
}
/**
For the output we evaluate the (obvious) total number of grid cells (*n*), the wall-clock time it took to perform 1000 timesteps (*perf.t*) and evaluate the error in the solution,
$$\eta=\int_{z=0}^{z=L0} \|v[\_]-v\| + \|u[\_]-u\| \mathrm{d}z$$. 
These three variables are written to a file named *ekmanfixed.dat*. 
*/
event eval(t=(10)){
  static FILE * fp = fopen("ekmanfixed.dat","w");
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
    fprintf(fpo,"%g\t%g\t%g\t%g\t%g\t%g\n", x,Delta,fabs(u[]-solu(x,Delta)),fabs(v[]-solv(x,Delta)),wv[],wv[]);
  fflush(fpo);
}

/**
We can view the results:

~~~gnuplot
  set xr [50:40000]
  set yr [0.000001:1]
  set logscale y
  set logscale x
  set xlabel 'Grid Cells'
  set ylabel 'Total Error'
    set size square
    plot    (500*x**(-2)) lw 3 lc rgb 'purple' title '{/Symbol \265}{N}^{-2}',\
            'ekmanfixed.dat' using  1:2 pt 4 title 'Total error'   
~~~

Concluding that the implementation of the solver and the boundary conditions are characterized by a second-order spatial accuracy. We also look at the error vs $\Delta$. 

~~~gnuplot This is the same data as Fig. 1a in the publication. 
  reset
  set logscale xy
  set xlabel 'Delta'
  set ylabel 'Error in u'
    set size square
    set key off
    plot  'ekmanfgoverview' u 2:3 pt 4    
~~~

*/