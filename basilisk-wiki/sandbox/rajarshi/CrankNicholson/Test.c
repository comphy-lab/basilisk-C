/**
# Poisson Solver on a Adaptive 1D Grid with Neumann Boundary Conditions and a Compact Polynomial Function

The second order Poisson scheme is used to solve a familiar compact function case with Neumann Boundary Conditions.
This case however uses adaptive grids which is done through the adapt_wavelet function.

*/

#include "grid/multigrid.h"
#include "utils.h"
#define BGHOST 1
#define Value 0

scalar A[],B[],Error[];

void solve(int depth){

   origin(-2,-2);
   L0 = 4;
   init_grid(1<<depth);

   foreach(){
      if(x*x<= 1 && y*y <= 1)
         A[] = (-pow(y,10)+5.*pow(y,8)-10.*pow(y,6)+10.*pow(y,4)-5.*pow(y,2)+1.)*(-pow(x,10)+5.*pow(x,8)-10.*pow(x,6)+10.*pow(x,4)-5.*pow(x,2)+1.);
      else
        A[] = 0.;
     
      B[] = 0.;
      Error[] = 0.;
     }

   foreach_dimension(){
      A[left] = neumann(Value);
      A[right] = neumann(Value)
      B[left] = neumann(Value);
      B[right] = neumann(Value);
     }

   boundary({A,B});
   
   foreach()
     B[] = (A[1,0]+A[0,1]-4*A[0]+A[-1,0]+A[0,-1])/(sq(Delta)); 

   foreach(){ 
     if(x*x<=1 && y*y<=1){
       Error[] = B[] - (-90.*pow(x,8) + 280.*pow(x,6) - 300.*pow(x,4) + 120.*pow(x,2) - 10.)*(-pow(y,10) + 5.*pow(y,8) - 10.*pow(y,6) + 10.*pow(y,4) - 5.*pow(y,2) + 1.) - (-90.*pow(y,8) + 280.*pow(y,6) - 300.*pow(y,4) + 120.*pow(y,2) - 10.)*(-pow(x,10) + 5.*pow(x,8) - 10.*pow(x,6) + 10.*pow(x,4) - 5.*pow(x,2) + 1.); 
      }
     else
       Error[] = B[];
    }

   FILE *fp1 = fopen("Error.dat","w");
   FILE *fp2 = fopen("Laplacian.dat","w");
   FILE *fp3 = fopen("ErrorvsGrid.dat","a");
   foreach(){
     fprintf(fp1,"%g %g %g \n",x,y,Error[]);
     fprintf(fp2,"%g %g %g \n",x,y,B[]);
   }
   fprintf(fp3,"%g %g \n",log10(pow(2,depth)),log10(normf(Error).max));
   fclose(fp1);
   fclose(fp2);
   fclose(fp3);
   printf("\n\n The Max Error is %g ",normf(Error).max);
}

int main(){
  system("rm -f ErrorvsGrid.dat");
  for(int depth=6; depth<=10;depth++)
     solve(depth);
}

/**
## Results
~~~gnuplot Error Distribution of the 4th order Laplacian scheme on a homogenous grid

set output 'ErrorDistribution.png'
splot 'Error.dat'
~~~

~~~gnuplot  of the Poisson 4th order solver on a heterogenous grid

set output 'Laplacian.png'
splot 'Laplacian.dat'
~~~

~~~gnuplot Error Scaling of the Poisson 4th order solver
set output 'ErrorScaling.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 1:2 via a, b
title_f(a,b) = sprintf('E(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 1:2 w p t 'ErrorScaling', f(x) t title_f(a,b)
~~~
*/