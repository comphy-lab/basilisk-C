/**
# Diffusion of a Gaussian Pulse on a locally-refined Grid
As a follow-up of [this](diffgausspulse.c) page, we study the effect of local grid refinement on the accuracy of the solution. On the aforementioned page we found that increasing the spatial resolution from 64 to 2048 grid point made the result much more accurate. This effect is most appreciable in the region of the Gaussian pulse core. Therefore we now run the simulation with a locally-refined grid.
*/
#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

double f(double t,double to,double xi,double D,double Delta){
  return (pow((to)/(to+t),0.5)*pow(M_PI*D*(t+to),0.5)*(erf((xi+(Delta/2))/(2*pow(D*(t+to),0.5)))-erf((xi-(Delta/2))/(2*pow(D*(t+to),0.5)))))/Delta;
}
double t0=1./3.;
double D=0.01;
const face vector DD[]= {D};
scalar c[];
int m,mm;
FILE * fp1;
char name1[100];
/**
Many different simulations are run. The grid is varied from $64/128$ to $1024/2048$ combined resolutions. We run with tree different settings.

1. Default Basilisk settings
2. Injection prolongaiton
3. Linear prolongation

*/
int main(){
  L0=5;
  X0=-L0/2;
  for (mm=0;mm<3;mm++){
    sprintf(name1,"toterrvsD%d.dat",mm);
    if (mm>0)
      fclose(fp1);
    fp1=fopen(name1,"w");
    for (m=0;m<5;m++)
      run(); 
  }
}

event init(t=0){
  DT=1e-4;
  init_grid(1<<(6+m)); // Initialize coarse grid
  refine(fabs(x)<0.4 && level<(7+m)); // Refine in the centre of the domain by one level
  if (mm==1)
    c.prolongation=refine_injection;
  if (mm==2)
    c.prolongation=refine_linear;
  foreach()
    c[]=f(t,t0,x,D,Delta);
  boundary({c});
}

event diffn(i++;t<=1.) {
  dt=dtnext(DT);
  diffusion(c,dt,DD); 
}

event end(t=1){
  char name[100];
  sprintf(name,"%dGauss%d.dat",mm,m);
  FILE * fp = fopen(name,"w");
  foreach()
    fprintf(fp,"%g\t%g\t%g\t%g\t%g\n",x,Delta,c[],f(t,t0,x,D,Delta),fabs(c[]-f(t,t0,x,D,Delta)));
  fflush(fp);
  double toterr=0; 
  foreach()
    toterr+=fabs(c[]-f(t,t0,x,D,Delta))*Delta;
  fprintf(fp1,"%d\t%g\t%g\n",m,L0/((double) pow(2.,(double) (m+6.))),toterr);
  fflush(fp1);
}
/**
## Results
Lets check if second-order-accurate scaling is enherited by locally refined grids:
~~~gnuplot
  set xr [0.002 : 0.1]
  set yr [0.000002 : 0.01]
  set logscale y 
  set logscale x
  set xlabel '{/Symbol D}_{max}'
  set ylabel 'Total Error'
  set key left top box 3
  set size square
    plot    (0.15*x**2) lw 3 lc rgb 'purple' title '{/Symbol \265}{/Symbol D}^{2}' ,\
            (0.05*x**1.2) lw 3 lc rgb 'green' title '{/Symbol \265}{/Symbol D}^{1.2}' ,\
            'toterrvsD0.dat' using 2:3 pt 4 title 'Default settings',\
            'toterrvsD1.dat' using 2:3 pt 4 title 'Injection prolongation',\
            'toterrvsD2.dat' using 2:3 pt 4 title 'Linear prolongation'
         ~~~
         
Even for these anisotropic grids, the second order scaling is maintianed when using the default settings. However, the injection method for prolongation does introduce an error that that does not scale second-order accurate. The convergence rate for the injection prolongation runs reflects a mix of first-order accurate and second order accurate formulations. Furthermore, A keen eye might spot that the linear prolongation is somewhat better than the default method. However, it does not give rise to a higher-order-convergence rate.  

## Analysis: Why is $refine\_linear$ more accurate?
When the coarse level solution is prolongated to the fine-level-solution halos at the resolution boundaries, the $c.prolongation$ attribute is used. In general, it is beneficial that numerical formulations enherit certain properties of the studied physical system. The diffusion equation we are solving is conservative for example. However, the default $refine\_bilinear$ technique is not, whearas the $refine\_linear$ is conservative. This does not mean that this technique is therefore by definition more accurate (e.g. $refine\_injection$ is also conservative but performs very poor for prolongation purposes). 

In order to answer the question of this subsection, let us approach the problem for another angle. Let us assume that the solution can be locally approximated by a second-order Taylor expansion arround $x_i$, 

$$c(x)\approx c_i+b(x-x_i)+a(x-x_i)^2+O\left((x-xi)^3\right),$$

where we have used $c(x_i)=c_i$ as notation. We choose to base our prolongation on this taylor expansion. This entails finding coeficients for $a$ and $b$ (and $c_i$?). We can solve a linear system for this by demanding that,

$$c_n= \frac{1}{\Delta} \int_{x_n-\Delta/2}^{x_n+\Delta/2} c(x)dx,$$

for $n=i+\left[-1,0,1\right]$. With $a,b$ and $c_i$ we can find conservative higher level solution ( $c^{\mathrm{high}}_j$ ) based on our Taylor expansion by definening it as,

$$c^{\mathrm{high}}_j=\frac{1}{\Delta^{high}}\int_{x^{\mathrm{high}}_j-\Delta^{high}}^{x^{\mathrm{high}}_j-\Delta^{high}}c(x)dx$$

If one does the math, it will turn out to be simply applying weights to the local 3-cell coarse level stencil. And these weights will turn out to be: 

~~~terminal
foreach_child()
    s[]=8.*(coarse(s,0)+coarse(s,child.x)-coarse(s,-child.x))/8;
~~~

This is exactly the same as what $refine\_linear$ does, but written down differently. Note that this similarity will not be inherited by higher dimensions. Therefore, I conclude that by coincidence (i.e. due to the limited degrees of freedom) $refine\_linear$ is third order accurate in one dimension. Therefore, it outperformed the second and first-order accurate (i.e. $bilinear$ and $injection$, respectively) formulations. Because we did not change the number of resolution boundaries in the runs on this page, we do not observe any different scaling compared to the second-order accurate formulation.  

## The next step
This manually (static) refined grid has learned us something on the importance of the threatment of resolution boundaries. However, many like to employ solution-adaptive refinement of their grid. Hence, we will study this as a next step [here](diffgausspulse3.c).
*/