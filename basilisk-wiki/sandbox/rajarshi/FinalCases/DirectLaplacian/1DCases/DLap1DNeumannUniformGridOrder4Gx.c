/**
## 1D Test case of 4th Order Laplacian operator on a Uniform Grid

In this code the error scaling of the Laplacian operator has been demonstrated on a Uniform Grid. A compact function which includes a polynomial function of order 10 is used.
*/

#include "grid/multigrid.h"
#include "utils.h"
#define BGHOST 2
#define Value 0
#define dimension 1

scalar A[],B[],Error[];

void solve(int depth){

   origin(-2);
   L0 = 4;
   init_grid(1<<depth);

   foreach(){
      if(x*x<=1){
        A[] = -pow(x,10) + 5.*pow(x,8) - 10.*pow(x,6) + 10.*pow(x,4) -5.*pow(x,2) + 1;
        }
      else
        A[] = 0.;
     
      B[] = 0.;
      Error[] = 0.;
     }

   foreach_dimension(){
      A[left] = neumann(Value);
      B[left] = neumann(Value);
     }

   boundary({A});
   
#if TREE
   face vector g[];
   foreach_face()
      g.x[] = (a[-2] - 27*a[-1] + 27*a[0] - a[1])/(24*Delta);
   boundary_flux({g});
   foreach()
     B[] = (g.x[1]-g.x[])/(Delta);

#else

   foreach(){
       B[] = (-A[-2]+16.*A[-1]-30.*A[0]+16.*A[1]-A[2])/(12*sq(Delta));
    }
#endif

   foreach(){ 
     if(x*x<=1){
       Error[] = B[] + 90.*pow(x,8) - 280.*pow(x,6) + 300.*pow(x,4) - 120.*pow(x,2) + 10.;
      }
      else
       Error[] = B[];
    }

   FILE *fp1 = fopen("Error.dat","w");
   FILE *fp2 = fopen("Laplacian.dat","w");
   FILE *fp3 = fopen("ErrorvsGrid.dat","a");
   foreach(){
     fprintf(fp1,"%g %g \n",x,Error[]);
     fprintf(fp2,"%g %g \n",x,B[]);
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
set xlabel 'x'
set ylabel 'Error'
set grid
plot 'Error.dat' u 1:2 w p t 'Error Lap(g(x))'
~~~

~~~gnuplot  Computed Laplacian

set output 'Laplacian.png'
set xlabel 'x'
set ylabel 'Error'
set grid
plot 'Laplacian.dat' u 1:2 w p t 'Laplacian'
~~~

~~~gnuplot Error Scaling of the Laplacian 4th order scheme
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