/**
# The Navier-Stokes solver is second-order accurate in time
Let us test that for a 1D viscous-flow problem, evaluated on a 2D grid.

##set-up
The set-up consists of a periodic (left-right) channel with no-slip boundary conditions (top-bottom). The initialized flow is:

$$ \overrightarrow{u_0}(x,y)= (u_0,v_0) =(\mathrm{cos}(y),0),$$

for $\mathrm{abs}(y) \leq \pi/2  $. Due to the fluid's viscousity ($\nu$), momemtum will diffuse over time, according to the (Navier-)Stokes equation, that is solved by;

$$\overrightarrow{u}(t)= \overrightarrow{u_0}e^{-\nu t}.$$

This equation will help to determine the error due to the time integration. 
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

const face vector muv[]={1.,1.};
double timestepp;
int j;

int main(){
  /**
  The grid resolution is chosen so fine that the presented errors here are dominated by the time stepping and are not significantly affected by the spatial discretization.
  */
  init_grid(1<<9);
  periodic(left);
  u.t[top]=dirichlet(0);
  u.t[bottom]=dirichlet(0.);
  mu=muv;
  L0=M_PI;
  X0=Y0=-L0/2.;
  /**
  Six experiments are run, the zeroth one ($j=0$) is not used because the timestepper thinks it has made a timestep with $\Delta t = 0$ for the -1-th timestep. 
  */
  for (j=0;j<7;j++){
    timestepp=2.*pow(0.5,(double)j);
    run();
  }
}

event init(t=0){
  CFL=10000.;
  foreach()
    u.x[]=cos(y);
  boundary(all);
  DT=timestepp;
}
/**
This event is used to let the run have an equidistant timestep for $j>0$.
*/
event setDT(i++){
if (i==0){
    DT=1.1*timestepp/0.1;
  }else{
    DT=timestepp;
  }
}
/**
## Output
After a single '$1/e$' timescale, the error in the numerically obtained solution is evaluated.
*/
event check(t=1.){
  static FILE * fp = fopen("resultvisc.dat","w");
  double err=0.;
  foreach()
    err+=fabs(u.x[]-cos(y)*exp(-t))*sq(Delta);
  if (j>0)
    fprintf(fp,"%d\t%g\n",i,err);
}  

/**
## Results

~~~gnuplot A first-order accurate (global) timestepping can be diagnosed? 
  set xr [0.5:40]
  set yr [0.02:1]
  set logscale y
  set logscale x
  set xlabel '{Used number of Timesteps}'
  set ylabel 'Total Error'
  set key box 1
  set size square
    plot    (1*x**-1) lw 3 lc rgb 'purple' title '∝{dt}^{-1}',\
            'resultvisc.dat' using 1:2 pt 4 title 'Error at t=1' 
~~~


*/

