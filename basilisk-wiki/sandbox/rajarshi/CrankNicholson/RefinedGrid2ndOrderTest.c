/**
# Poisson 2nd Order Solver - Error Scaling

We solve the Poisson equation using a 2nd order Finite Volume Solver.
This exercise is performed to test the Error scalability of the solver for which
we use differently resolved grids & study how the maximum error behaves */


#include "grid/tree.h"
#include "poisson.h"
#include "utils.h"
#define BGHOSTS 2

void Solve(int depth){
  
   clock_t start,end;
   start = clock();
   L0=1;
   init_grid(1);
   refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));
   periodic(left);
   periodic(top);

   scalar A[], B[], Res[], dA[];
   foreach(){
     //B[] = -2.*(sin(2.*pi*(x+Delta/2.))-sin(2.*pi*(x-Delta/2.)))*(sin(2.*pi*(y+Delta/2.))-sin(2.*pi*(y-Delta/2.)))/sq(Delta);
       B[] = -8.*pi*pi*cos(2.*pi*x)*cos(2.*pi*y);
       A[] = 0.;
    }
   boundary({A,B});
   double maxres = 1E-04;
   mgstats Sol = poisson(A,B, tolerance = maxres);
   stats s = statsf(A);
   foreach()
     A[] = A[] - s.sum/s.volume;
   scalar Error[];
   foreach()
      Error[] = A[] - cos(2.*pi*x)*cos(2.*pi*y);
   static FILE *fp,*fp2;
   fp = fopen("ErrorvsGrid.dat","a");
   fp2 = fopen("Time2nd.dat","a");
   fprintf(fp,"%d %g %g %g %g %g \n",depth,pow(2,depth),normf(Error).max,log10(pow(2,depth)),log10(normf(Error).max),s.sum/s.volume);
   fclose(fp);

   printf("\n The Depth is = %d \n\n",depth);
   printf("\n No. of Iterations = %d  and the Residual = %.8f \n\n",Sol.i,Sol.resa);
   
   printf("\n The Error is = %g \n\n\n",normf(Error).max);
   end = clock();
   printf("The total time is : %g \n\n\n\n\n\n",(end-start)/(double)CLOCKS_PER_SEC);
   fprintf(fp2,"%g    %g     %d      %g\n",pow(2,depth),(end-start)/(double)CLOCKS_PER_SEC,Sol.i,Sol.resa);
   fclose(fp2);
}

int main(){
 remove("ErrorvsGrid.dat");
 remove("Time2nd.dat");
 int depth;
 for(depth=6;depth<=8;depth++)
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
fit f(x) 'ErrorvsGrid.dat' u 4:5 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 4:5 w p t 'ErrorScaling', f(x) t title_f(a,b)
~~~

~~~gnuplot Iterations of the Poisson 4th order solver
set output 'Iterations.png'
set xlabel 'GridPoints'
set ylabel 'No. of Iterations'
plot 'Time2nd.dat' u 1:3 w p t 'Iterations'
~~~

*/