/**
# Convergence of the Poisson solver using an adaptive 1D grid
Very similar to the spatial-convergence test of the Poisson solver on a [regular grid](testpoisson.c), we now look at the convergence properties of the solver with the refinement criterion, $\zeta$. For this case we adapt our grid based on the estimated error in the source term. 

## Method
We check how the error scales with the refinement criterion. We limit the maximum resolution to an $N=4096$ grid. Hence convergence is only expected with respect to the solution obtained from using a regular grid, using that numer of grid cells.  

## Code:
*/
#include "grid/bitree.h"
#include "poisson.h"
#include "utils.h"

#define sol(x,c1,c2) ((c1*x)+c2+(pow(M_PI,0.5)*(x)*erf(x)/2)+(exp(-(x*x))/2))

double c1,c2,err;
scalar b[],a[];
mgstats mg;
FILE * fp1;
FILE * fp2;
char name[100];
double errzero;

double numintsol(double tol,double xi, double D,double c1,double c2){
  double into=5;
  double intn=10;
  double integral;
  double j=10;
  while ((fabs(into-intn)/fabs(intn))>tol){
    into=intn;
    intn=0.;
    integral=0;
    for (int m=0;m<j;m++){
      double xp=(xi-(D/2)+(D/(j*2))+(((double)m)*D/j));
      integral+=sol(xp,c1,c2);
    }
    intn=integral/j;
    j=j*2;
  }
  return intn;
}

double f(double xi,double Delta){
  return (pow(M_PI,0.5)*(erf(((xi)+(Delta/2)))-erf(((xi)-(Delta/2)))))/(2*Delta);
}

int main(){
  L0=10;
  X0=-L0/2;
  a[right]=dirichlet(0.);
  a[left]=neumann(0.);
  c1 = pow(M_PI,0.5)/2;
  c2 =-sol(5,c1,0);
  fp2=fopen("Convergencewithzeta","w");
  TOLERANCE=10; //Set large Tolerance
  int n=0;
  /**
  We now have a loop where we iteratively set a smaller refinement criterion ($\zeta$)
  */
  for (int j=0;j<20;j++){ 
    double zeta=1./pow(2.,(double)j);
    if (j==0) //Benchmark result
      zeta=0.;
    init_grid(1<<12);
    foreach(){
      a[]=0.;
      b[]=f(x,Delta);
    }
    while (adapt_wavelet({b},(double[]){zeta},12).nc){
      foreach(){
	a[]=0.;
	b[]=f(x,Delta);
      }
      boundary({a,b});
    }
    boundary(all);
    sprintf(name,"MGcycles%d.dat",j);
    fp1=fopen(name,"w");
    fprintf(ferr,"Grid j=%d\n",j);
    boundary({a,b});
    /**
    Ten Multigrid cycles of the solver are performed. We learned [here](testpoisson.c) that this should be enough to reach the converged solution for each grid.  
    */
    for (int i=0;i<10;i++){
      mg = poisson(a,b); //Solve the system
      err=0;
      double errz=0;
      n=0;
      foreach(){
	err+=fabs(a[]-numintsol(10e-5*pow(Delta,3.),x,Delta,c1,c2))*Delta;
	errz+=fabs(a[]-sol(x,c1,c2))*Delta;
	n++;
      }
      if(j==0)
	errzero=err;
      fprintf(fp1,"%d\t%d\t%g\t%g\t%g\t%d\t%g\t%g\n",i,mg.i,mg.resb,mg.resa,mg.sum,mg.nrelax,err,errz);
    }
    fclose(fp1);
    if(j>0)
      fprintf(fp2,"%d\t%g\t%g\t%d\n",j,zeta,fabs(err-errzero),n);
    fflush(fp2);
  }
  fclose(fp2);
}
/**
## Results
First we check the correlation between the number of grid cells and the refinement criterion. 

~~~gnuplot Number of grid cells versus the refinement criterion
set xr [0.000001:1]
set yr [1:4000]
set logscale y
set logscale x
set xlabel '{/Symbol \z}'
set ylabel 'Number of Grid Cells'
set key off
plot 'Convergencewithzeta' using 2:4 pt 3 
~~~

As expected, a lower refinment criterion corresponds to a more refined grid. 

~~~gnuplot Convergence of the solver with the refinement criterion 

set yr [0.0000001:1]

set ylabel '{Error}'

plot 'Convergencewithzeta' using 2:3 pt 4
~~~

The error convergence does not seems to be very well behaved for $\zeta>0.001$. 

Additionally, we check how the residual and the error correlate during the ten succesive multigrid cycles we have performed on each grid (i.e for each $\zeta$).
~~~gnuplot
  
  set xr [1e-13:1]
  set yr [1e-8:10]
  set logscale x
  
  set xlabel 'Residual'
  set ylabel 'Error'
  set key left top box 1
  set size square
  plot       'MGcycles8.dat' using 4:8 with lines linecolor rgb "#0000ff" lw 3 title '{/Symbol \z}=0.25' ,\
     'MGcycles2.dat' using 4:8 with lines linecolor rgb "#0033bb"  lw 3 title '{/Symbol \z}=1/16' ,\
  'MGcycles4.dat' using 4:8 with lines linecolor rgb "#0066bb"  lw 3 title '{/Symbol \z}=1/64' ,\
  'MGcycles6.dat' using 4:8 with lines linecolor rgb "#006699"  lw 3  title '{/Symbol \z}=1/256' ,\
  'MGcycles8.dat' using 4:8  with lines linecolor rgb "#009999"  lw 3 title '{/Symbol \z}=1/1024' ,\
  'MGcycles10.dat' using 4:8 with lines linecolor rgb "#009966"  lw 3  title '{/Symbol \z}=1/4096' ,\
  'MGcycles12.dat' using 4:8 with lines linecolor rgb "#00bb66"  lw 3  title '{/Symbol \z}=1/16384' ,\
   'MGcycles14.dat' using 4:8 with lines linecolor rgb "#00bb33"  lw 3  title '{/Symbol \z}=1/64436' ,\
  'MGcycles16.dat' using 4:8 with lines linecolor rgb "#00dd33"  lw 3 title '{/Symbol \z}=1/262144' ,\
  'MGcycles18.dat' using 4:8 with lines linecolor rgb "#00ff00"  lw 3 title '{/Symbol \z}=1/1048576' ,\
  'MGcycles0.dat' using 4:8 with lines linecolor rgb "#dd33dd" lw 4 title 'N=4096 grid'
 
~~~

I am not sure how to interpret this plot. But it seems to contain some interesting information. It atleast shows that the choice of using ten Multigrid cycles was indeed enough for all grids.
*/