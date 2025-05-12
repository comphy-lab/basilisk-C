/**
# Poisson 4th Order Solver - Error Scaling

We solve the Poisson equation using a 4th order Finite Volume Solver.
This exercise is performed to test the Error scalability of the solver on a periodic domain for which
we use differently resolved grids & study how the maximum error behaves

*/

// #include "grid/multigrid.h"
#if 1
# include "Poisson-helmholtz4th.h"
#else
# include "poisson.h"
#endif
#include "utils.h"
#define BGHOSTS 2

void Solve(int depth){
  
   clock_t start,end;
   start = clock();
   L0=1;
   init_grid(1<<depth);
   periodic(left);
   periodic(top);

   scalar A[], B[], Res[], dA[];

  /**
   A known solution of the poisson problem is chosen so the numerical solution can be later compared 
   with the analytical solution
*/
  
  foreach(){
     B[] = -2.*(sin(2.*pi*(x+Delta/2.))-sin(2.*pi*(x-Delta/2.)))*(sin(2.*pi*(y+Delta/2.))-sin(2.*pi*(y-Delta/2.)))/sq(Delta);
     A[] = 0.;
    }
   boundary({A,B});
  
/**
   The maximum residual of 1E-06 is chosen and the calls to the fourth order poisson solver is made
   The header file for the same is saved in Poisson-helmholtz4th.h
*/
  
   double maxres = 1E-09;
   mgstats Sol = poisson(A,B, tolerance = maxres);
   scalar Error[];
   foreach()
      Error[] = A[] - cos(2.*pi*x)*cos(2.*pi*y);

   fprintf(stderr,"%d %g %g %g %g \n",depth,pow(2,depth),normf(Error).max,log10(pow(2,depth)),log10(normf(Error).max));

   printf("\n The Depth is = %d \n\n",depth);
   printf("\n No. of Iterations = %d  and the Residual = %.8f \n\n",Sol.i,Sol.resa);
   
   printf("\n The Error is = %g \n\n\n",normf(Error).max);
   end = clock();
   printf("The total time is : %g \n\n\n\n\n\n",(end-start)/(double)CLOCKS_PER_SEC);

}

/**
  From the Main function successive calls are made to the solve function with different Grid densities
  The errors are stored in a File and later plotted to show a slope variation
*/

int main(){
  
 int depth;
 for(depth=3;depth<=9;depth++)
    Solve(depth);  
}

/**
## Results
~~~gnuplot Error Scaling of the Poisson 4th order solver
set output 'ErrorScaling.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'log' u 4:5 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'log' u 4:5 w p t 'ErrorScaling', f(x) t title_f(a,b)
~~~


*/
