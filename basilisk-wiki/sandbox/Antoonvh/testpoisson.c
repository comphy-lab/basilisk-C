/**
# Convergence of the Poisson solver in 1D
On this page we try to find some of the convergence properties of the Basilisk poisson solver. The used solver is of the 'multi-grid' type and untilizes an iterative scheme to arrive at an converged solution. By default the solver iterates untill the maximum residual is below some set tolerance. On this page we solve a Poisson problem on 11 different grids with various levels of refinement and let the solver iterate (i.e. Cycle) ten times. We will check the convergence properties with respect to the cycles applied and the resolution. The goal is to check if we obtain sensible results so we can extend our analysis to grids generated according to an adaptive algorithm. 

## The chosen Problem
The Poisson problem we will be solving for reads,

$$\frac{\partial^2 a}{\partial x^2}=e^{-x^2}.$$

According to [Wolfram|Alpha](https://www.wolframalpha.com/input/?i=int+int+exp(-x%5E2)) the corresponding solution is,

$$a=c_2+c_1x+\frac{\sqrt{\pi}}{2}\mathrm{erf}(x)+\frac{e^{-x^2}}{2},$$

with constants $c_1$ and $c_2$ determined by the boundary conditions. 

## The script
The Poisson solver is included and we opt for a one dimensional tree grid. The tree functionality will be usefull later. 
*/
#include "grid/bitree.h"
#include "poisson.h"
/**
The analyical solution can be evaluated to check the numerical solution. Furthermore, we initialize some usefull stuff. 
*/
#define sol(x,c1,c2) ((c1*x)+c2+(pow(M_PI,0.5)*(x)*erf(x)/2)+(exp(-(x*x))/2))

double c1,c2;
scalar b[],a[];
mgstats mg;
FILE * fp1;
FILE * fp2;
char name[100];
/**
Since our solution represents volume averaged values we need to translate the analytical solution to a locally averaged one. Unfortunately, Wolfram Alpha was unable to provide me with the the integal form of the solution. Therefore a numerical integrator of the analytical solution is defined
*/
double numintsol(double tol,double xi, double D,double c1,double c2){ //A Riemann integrator
  double into=5;
  double intn=1;
  double integral;
  double j=10;
  while ((fabs(into-intn)/fabs(intn))>tol){ // Perform the integration untill the integral has converged
    into=intn;
    intn=0.;
    integral=0;
    for (int m=0;m<j;m++){ // Summate the analytical solution j times 
      double xp=(xi-(D/2)+(D/(j*2))+(((double)m)*D/j)); //At equally spaced locations in the grid box
      integral+=sol(xp,c1,c2);
    }
    intn=integral/j; 
    j=j*2; // Increase j if the integral has not converged upto the set standards
  }
  return intn;
}
/**
Similarly, the source term should be defined consistently.
*/
double f(double xi,double Delta){
  return (pow(M_PI,0.5)*(erf(((xi)+(Delta/2)))-erf(((xi)-(Delta/2)))))/(2*Delta);
}
/**
The rest of the script occurs in the $main()$ function. The spatial extend of the grid is defined, suitable boundary conditions are chosen and the corresponding values for $c_1$ and $c_2$ are set. 
*/
int main(){
  L0=10;
  X0=-L0/2;
  a[right]=dirichlet(0.);
  a[left]=neumann(0.); //Default Basilisk
  c1=pow(M_PI,0.5)/2;
  c2 =-sol(5,c1,0);
  fp2=fopen("gridcycles.dat","w");

  /**
  A loop is used to iterate over 11 different resolutions, varying from $N=8$  to $N=8192$ grid cells. Each time the solution and source term are initialized.
  */
  for (int j=0;j<11;j++){
    init_grid(1<<(j+3));
    sprintf(name,"MGcycles%d.dat",N);
    fp1=fopen(name,"w");
    foreach(){
      a[]=0.; // This choice is the problem of the poisson solver 
      b[]=f(x,Delta);
    }
    double err;
    TOLERANCE=10; // Large Tolerance to prevent unwanted iterations
    fprintf(ferr,"Grid N=%d\n",N);
    boundary({a,b}); 
    /**
    On each grid we let the Poisson solver iterate 10 times. After every iteration we write down the statistics of the solver and the error in the obtained solution.
    */
    for (int i=0;i<10;i++){
      mg = poisson(a,b); //Solve the system
      err=0;
      double err2=0;
      foreach(){
        err+=fabs(a[]-numintsol(10e-4*pow(Delta,3.),x,Delta,c1,c2))*Delta; //We evaluate the analytical solution with 3-rd order accuracy
        err2+=fabs(a[]-sol(x,c1,c2))*Delta; //We do not use this
      }
      fprintf(fp1,"%d\t%d\t%g\t%g\t%g\t%d\t%g\t%g\n",i,mg.i,mg.resb,mg.resa,mg.sum,mg.nrelax,err,err2);
    }
    fclose(fp1);
    fprintf(fp2,"%d\t%g\t%g\n",j,L0/((double)N),err);
    fflush(fp2);
  }
  fclose(fp2);
}

/**
## Results
First we check how the maximum residual decreases with every iteration for three grids with different resolutions.

~~~gnuplot
  
  set xr [-0.5:9.5]
  set yr [1e-15:2]
  set logscale y
  
  set xlabel 'Iteration'
  set ylabel 'Residual'
  set key right top box 1
  set size square
  plot     (5**(-2*x)) lw 8 lc rgb 'purple' title '{/Symbol \265}c^{-2i}' ,\
  'MGcycles32.dat' using 1:3 pt 4 title ' N=32' ,\
  'MGcycles512.dat' using 1:3 pt 3 title ' N=512',\
  'MGcycles8192.dat' using 1:3 pt 3 title ' N=8192'

~~~

The convergence regions seems to start from $i=0$ and end at a iteration value that is resolution dependent. For the higher resolutions the residual seems to stop converging earlier and at a with a higher residual than the coarse grid runs. 

Next we check how the error in the solution,

~~~gnuplot
  
  set xr [-0.5:9.5]
  set yr [1e-8:2]
  set logscale y
  
  set xlabel 'Iteration'
  set ylabel 'Error'
  set key right top box 1
  set size square
  plot     (5**(-2*x)) lw 8 lc rgb 'purple' title '{/Symbol \265}c^{-2i}' ,\
  'MGcycles32.dat' using 1:8 pt 4 title ' N=32' ,\
  'MGcycles512.dat' using 1:8 pt 3 title ' N=512',\
  'MGcycles8192.dat' using 1:8 pt 3 title ' N=8192'
~~~

Now we see that the smaller residuals of the coarse runs do not translate into smaller errors in the solution. Also we see that with the first few iterations the fine grid solutions are not more accurate. It seems that only the fine-grid runs benefit from the large number of iterations. However, even for the finest-grid run the solution converges before the 10-th iteration. from We can check this  converged solution error in a bit more detail.

~~~gnuplot
  
  set xr [0.001:2]
  set yr [1e-8:1]
  set logscale x
  
  set xlabel '{/Symbol D}'
  set ylabel 'Error'
  set key left top box 1
  set size square
  plot     (0.1*(x)**2) lw 8 lc rgb 'purple' title '{/Symbol \265}{/Symbol D}^{2}' ,\
  'gridcycles.dat' using 2:3 pt 5 title 'Total error'
 
~~~

So we can conclude that the Poisson solver is second-order accurate (i.e for the converged solution). Notice that this convergence takes more iterations for the higher resolution runs.

For some elucidation, here is an additional plot.

~~~gnuplot
  
  set xr [1e-13:1]
  set yr [1e-8:10]
  set logscale x
  
  set xlabel 'Residual'
  set ylabel 'Error'
  set key right bottom box 1
  set size square
  plot       'MGcycles8.dat' using 4:8 with lines linecolor rgb "#0000ff" lw 3 title '      N=8' ,\
  'MGcycles16.dat' using 4:8 with lines linecolor rgb "#0033dd" lw 3 title 'N=16' ,\
  'MGcycles32.dat' using 4:8 with lines linecolor rgb "#0033bb"  lw 3 title 'N=32' ,\
  'MGcycles64.dat' using 4:8 with lines linecolor rgb "#0066bb"  lw 3 title 'N=64' ,\
  'MGcycles128.dat' using 4:8 with lines linecolor rgb "#006699"  lw 3  title 'N=128' ,\
  'MGcycles256.dat' using 4:8  with lines linecolor rgb "#009999"  lw 3 title 'N=256' ,\
  'MGcycles512.dat' using 4:8 with lines linecolor rgb "#009966"  lw 3  title 'N=512' ,\
  'MGcycles1024.dat' using 4:8 with lines linecolor rgb "#00bb66"  lw 3  title 'N=1024' ,\
   'MGcycles2048.dat' using 4:8 with lines linecolor rgb "#00bb33"  lw 3  title 'N=2048' ,\
  'MGcycles4096.dat' using 4:8 with lines linecolor rgb "#00dd33"  lw 3 title 'N=4096' ,\
  'MGcycles8192.dat' using 4:8 with lines linecolor rgb "#00ff00"  lw 3 title 'N=8192' 
 
~~~
*/