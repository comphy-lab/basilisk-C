/**
# Poisson Solver on a Uniform 1D Grid with Neumann Boundary Conditions and a Compact Polynomial Function

We solve the Poisson equation using and Order 2 Scheme.
*/

#include "grid/multigrid.h"
#include "poisson.h"
#include "Laplacian.h"
#include "utils.h"
#define BGHOSTS 1
#define dimension 1
#define Value 0

scalar A[], B[], Res[], dA[];

void Compute_B(){

 foreach(){
     if(x*x<=1)
        B[] = (-90.*pow(x,8) + 280.*pow(x,6) - 300.*pow(x,4) + 120.*pow(x,2) - 10.);
     else
        B[] = 0.;
     }
  
  stats s1 = statsf(B);
  foreach()
      B[] -= s1.sum/s1.volume;
  
  B[left] = neumann(Value);
  B[right] = neumann(Value);
  B.prolongation = refine_biquartic;
  boundary({B});

}

void Solve(int depth){
  
   L0=4;
   origin(-2);
   init_grid(1<<depth);
   
   Compute_B();
   foreach()   
     A[] = 0.;
    
   A[left] = neumann(Value);
   A[right] = neumann(Value); 
   boundary({A});

#define NITER 30
   int nrelax = 4;
   const scalar lambda[] = 0.;
   struct Poisson p;
   p.a = A; p.b = B; p.alpha = unityf; p.lambda = lambda;
   scalar * lres = {Res};

  for (int i = 0; i < NITER; i++) {
    residual ({A}, {B}, lres, &p);
    mg_cycle ({A}, lres, {dA}, relax, &p, nrelax, 0, depth());
    double max = 0.;
    foreach()
      if (fabs(Res[]) > max)
	max = fabs(Res[]);
    printf("\n The max Residual after Iteration no %d is %g (Depth = %d)",i,max,depth);
  }

   scalar Error[];
   foreach(){
     if(x*x<=1)  
       Error[] = A[] - (-pow(x,10)+5.*pow(x,8)-10.*pow(x,6)+10.*pow(x,4)-5.*pow(x,2)+1.);
     else
       Error[] = A[];   
   }
   
   stats s = statsf(Error);
   foreach(){
     Error[] -= s.sum/s.volume;
     A[] -= s.sum/s.volume;
    }
   printf("\n Hello");


   FILE *fp,*fp1,*fp2;

   fp1 = fopen("Solution.dat","w");
   foreach()
    fprintf(fp1,"%g %g\n",x,A[]);
   fclose(fp1);

   fp2 = fopen("ErrorDistribution.dat","w");
   foreach()
    fprintf(fp2,"%g %g\n",x,Error[]);
   fclose(fp2);

   fp = fopen("ErrorvsGrid.dat","a");
   fprintf(fp,"%g %g \n",log10(pow(2,depth)),log10(normf(Error).max));
   fclose(fp);

}

int main(){
  
 system("rm -f ErrorvsGrid.dat");
 for(int depth=6;depth<=10;depth++)
    Solve(depth);  
}

/**
## Results

~~~gnuplot Error Distribution for Residual operator
set output 'ErrorDistribution.png'
plot 'ErrorDistribution.dat' u 1:2 w p t 'Error_Distribution'
~~~

~~~gnuplot Solution Distribution
set output 'Solution.png'
plot 'Solution.dat' u 1:2 w p t 'Solution_Distribution'
~~~

~~~gnuplot Error Scaling of the Poisson 4th order solver
set output 'ErrorScaling.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 1:2 w p t 'ErrorScaling', f(x) t title_f(a,b)
~~~

*/