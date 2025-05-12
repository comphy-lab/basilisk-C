/**
# Poisson Solver on a Adaptive 1D Grid with Neumann Boundary Conditions and a Compact Polynomial Function

The second order Poisson scheme is used to solve a familiar compact function case with Neumann Boundary Conditions.
This case however uses adaptive grids which is done through the adapt_wavelet function.

*/

#include "grid/bitree.h"
#include "poisson.h"
#include "Laplacian.h"
#include "utils.h"
#define BGHOSTS 1
#define dimension 1
#define Value 0

scalar A[], B[], Res[], dA[];

/**
This is a function, which will be called whenever a re-computation of the Rhs is required. Mostly this will be needed for initialization and after every adaptive refinement step.
*/

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

/**
This is the Solve function, which starts with setting the domain, and initializing an uniform grid with a given depth. Then a call is made to initialize the RHS.
Boundary conditions are applied on the function a (i.e. LHS, which needs to be computed) and a biquartic prolongation function is chosen from the header file Laplacian.h.
The poisson problem is set up, and utilizing the functions residual and mg_solve in the header file poisson.h, a set of iterations are run to achieve a convergent residual.
These calculations are all on the uniform Multigrid.

Once a converged solution is obtained on the Multigrid a further series of iterations is run on the grid which is progressively adapted using the adapt_wavelet function.
The adapt wavelet targets the Scalar field 'a[]' with a tolerance given by cmax (which is a value passed on to the Solve function from the main function). The minlevel is
set as the same as the depth of the multigrid, and the maxlevel is set as 10 (without this no adaptation is observed) . We will control the simulation only with the help of the cmax values. The final bracket set illustrates the set of variables which needs to be interpolated after the adaptive refinement has taken place. This includes just 'a[]'. The 'b[]' is updated using the call to Compute_B()

*/


void Solve(double cmax){
  
   int depth = 4;
   double max;
   L0=8;
   origin(-4);
   init_grid(1<<depth);  
   Compute_B();

   foreach()   
     A[] = 0.;
    
   A[left] = neumann(Value);
   A[right] = neumann(Value); 
   A.prolongation = refine_biquartic;
   boundary({A});

#define NITER 5
#define ANITER 30

   int nrelax = 4;
   const scalar lambda[] = 0.;
   struct Poisson p;
   p.a = A; p.b = B; p.alpha = unityf; p.lambda = lambda;
   scalar * lres = {Res};

/**
 Iterations on a Uniform grid
*/
  
  for (int i = 0; i < NITER; i++) {
    residual ({A}, {B}, lres, &p);
    mg_cycle ({A}, lres, {dA}, relax, &p, nrelax, 0, depth());
    max = 0.;
    foreach()
      if (fabs(Res[]) > max)
	max = fabs(Res[]);
    printf("\n The max Residual after Iteration no %d is %g (Depth = %d)",i,max,depth);
  }

/**
 Iterations on an adaptive grid
*/
  
  for (int i = 0; i < ANITER; i++) {
    residual ({A}, {B}, lres, &p);
    mg_cycle ({A}, lres, {dA}, relax, &p, nrelax, 0, depth());
    max = 0.;
    foreach()
      if (fabs(Res[]) > max)
	max = fabs(Res[]);
    astats s = adapt_wavelet ({A}, (double []){cmax},minlevel=depth,maxlevel=10,{A});
    Compute_B();
    printf("\n \n The Iteration Count : %d , The Max Residual : %g",i,max);
    printf("\n The number of refined grids : %d , The number of coarsened grids : %d , Total Grids : %d ",s.nf,s.nc,grid->tn); 
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
   fprintf(fp,"%g %g \n",log10(grid->tn),log10(normf(Error).max));
   fclose(fp);

}

int main(){
  
 system("rm -f ErrorvsGrid.dat");
 for(double cmax=1E-03;cmax>=1E-05;cmax=cmax/10.)
    Solve(cmax);  
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