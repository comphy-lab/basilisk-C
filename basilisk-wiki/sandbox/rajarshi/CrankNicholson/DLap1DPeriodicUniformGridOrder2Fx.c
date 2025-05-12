/**
## 1D Test case of 2nd Order Laplacian operator on a Uniform Grid with Periodic Boundary Conditions

In this code the error scaling of the Laplacian operator has been demonstrated. A compact function which includes a polynomial function of order 8 is used.
*/

#include "grid/multigrid.h"
#include "utils.h"
#include "Laplacian.h"
#define BGHOST 2
#define Value 0


scalar A[],B[],Error[];

void solve(int depth){

   origin(-2);
   L0 = 4;
   init_grid(1<<depth);
   
   foreach(){
     A[] = sin(2.*pi*x);
     B[] = 0.;
     Error[] = 0.;
   }
 
   A.prolongation = refine_biquintic;
   periodic(right);
   boundary({A});
   
#if TREE
   face vector g[];
   foreach_face()
      g.x[] = (A[-2] - 15.*A[-1] + 15.*A[] - A[1])/(12.*Delta);
   boundary_flux({g});
   foreach()
     B[] = (g.x[1] - g.x[])/(Delta);

#else

  foreach()
     B[] = (-A[-2]+16.*A[-1]-30.*A[0]+16.*A[1]-A[2])/(12.*sq(Delta));

#endif

   foreach()
      Error[] = B[] + 4.*pi*pi*sin(2.*pi*x);

   FILE *fp1 = fopen("Error.dat","w");
   FILE *fp2 = fopen("Laplacian.dat","w");
   FILE *fp3 = fopen("ErrorvsGrid.dat","a");
   foreach(){
     fprintf(fp1,"%g %g %g %g \n",x,A[],B[],Error[]);
     fprintf(fp2,"%g %g \n",x,B[]);
   }
   fprintf(fp3,"%g %g \n",log10(pow(2,depth)),log10(normf(Error).rms));
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
plot 'Error.dat' u 1:4 w p t 'Error'
~~~

~~~gnuplot  of the Poisson 4th order solver on a heterogenous grid

set output 'Laplacian.png'
set xlabel 'x'
set ylabel 'Error'
set grid
plot 'Laplacian.dat' u 1:2 w p t 'Laplacian'
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