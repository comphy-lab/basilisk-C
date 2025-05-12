/**
# Poisson Solver 2nd Order Scheme on an Adaptive Grid

The second order Poisson scheme is used to solve a familiar compact function case with Neumann Boundary Conditions.
This case however uses adaptive grids which is done through the adapt_wavelet function.

*/

#include "grid/quadtree.h"
#include "poisson.h"
#include "Laplacian.h"
#include "utils.h"
#define BGHOSTS 1
#define Value 0

scalar a[], b[], res[], dp[], Error[];

/**
This is a function, which will be called whenever a re-computation of the RHS is required. Mostly this will be needed for initialization and after every adaptive refinement step.
*/

void Compute_B(){

 foreach() {
    
    if(x*x<=1 && y*y <=1)
      b[] = (((-10.*pow(x+Delta/2.,9)+40.*pow(x+Delta/2.,7)-60.*pow(x+Delta/2.,5)+40.*pow(x+Delta/2.,3)-10.*(x+Delta/2.))-(-10.*pow(x-Delta/2.,9)+40.*pow(x-Delta/2.,7)-60.*pow(x-Delta/2.,5)+40.*pow(x-Delta/2.,3)-10.*(x-Delta/2.))) * ((-pow(y+Delta/2,11)/11. + 5.*pow(y+Delta/2.,9)/9. -10.*pow(y+Delta/2.,7)/7. + 10.*pow(y+Delta/2.,5)/5. - 5.*pow(y+Delta/2.,3)/3. + (y+Delta/2.)) - (-pow(y-Delta/2,11)/11. + 5.*pow(y-Delta/2.,9)/9. -10.*pow(y-Delta/2.,7)/7. + 10.*pow(y-Delta/2.,5)/5. - 5.*pow(y-Delta/2.,3)/3. + (y-Delta/2.))) + ((-10.*pow(y+Delta/2.,9)+40.*pow(y+Delta/2.,7)-60.*pow(y+Delta/2.,5)+40.*pow(y+Delta/2.,3)-10.*(y+Delta/2.))-(-10.*pow(y-Delta/2.,9)+40.*pow(y-Delta/2.,7)-60.*pow(y-Delta/2.,5)+40.*pow(y-Delta/2.,3)-10.*(y-Delta/2.))) * ((-pow(x+Delta/2,11)/11. + 5.*pow(x+Delta/2.,9)/9. -10.*pow(x+Delta/2.,7)/7. + 10.*pow(x+Delta/2.,5)/5. - 5.*pow(x+Delta/2.,3)/3. + (x+Delta/2.)) - (-pow(x-Delta/2,11)/11. + 5.*pow(x-Delta/2.,9)/9. -10.*pow(x-Delta/2.,7)/7. + 10.*pow(x-Delta/2.,5)/5. - 5.*pow(x-Delta/2.,3)/3. + (x-Delta/2.))))/sq(Delta);
     
    else
      b[] = 0.;
  }
    
  stats s = statsf(b);
  foreach()
    b[] -= s.sum/s.volume;  
 
  b[left] = neumann(Value);
  b[right] = neumann(Value);
  b[top] = neumann(Value);
  b[bottom] = neumann(Value);
  b.prolongation = refine_biquartic;
  boundary ({b});

}

/**
This is the Solve function, which starts with setting the domain, and initializing an uniform grid with a given depth. Then a call is made to initialize the RHS.
Boundary conditions are applied on the function a (i.e. LHS, which needs to be computed) and a biquartic prolongation function is chosen from the header file Laplacian.h.
The poisson problem is set up, and utilizing the functions residual and mg_solve in the header file poisson.h, a set of iterations are run to achieve a convergent residual.
These calculations are all on the uniform Multigrid.
*/

void solve (double cmax)
{ 
  int depth = 4;
  origin (-4,-4);
  int nrelax = 4;
  L0 = 8;
  init_grid(1<<depth);
  Compute_B();
  foreach()
   a[] = 0.;
  
  a[left] = neumann(Value);
  a[right] = neumann(Value);
  a[top] = neumann(Value);
  a[bottom] = neumann(Value);
  a.prolongation = refine_biquartic;
  boundary({a,dp});

  #define NITER 10
  #define ANITER 60

  const scalar lambda[] = 0.;
  struct Poisson p;
  p.a = a; p.b = b; p.alpha = unityf; p.lambda = lambda;
  scalar * lres = {res};
  
  double max;
  for (int i = 0; i < NITER; i++) {
    max=0.;
    residual ({a}, {b}, lres, &p);
    mg_cycle ({a}, lres, {dp}, relax, &p, nrelax, 0, depth());
    foreach()
      if (fabs(res[]) > max)
	max = fabs(res[]);
    printf("\n \n The Iteration Count : %d , The Max Residual : %g",i,max);
    }

/**
Once a converged solution is obtained on the Multigrid a further series of iterations is run on the grid which is progressively adapted using the adapt_wavelet function.
The adapt wavelet targets the Scalar field 'a[]' with a tolerance given by cmax (which is a value passed on to the Solve function from the main function). The minlevel is
set as the same as the depth of the multigrid, and the maxlevel is left open. We will control the simulation only with the help of the cmax values. The final bracket set
illustrates the set of variables which needs to be interpolated after the adaptive refinement has taken place. This includes just 'a[]'. The 'b[]' is updated using the call to Compute_B()
*/
  
  for(int i=0; i<ANITER; i++){
    max=0.;
    residual ({a}, {b}, lres, &p);
    mg_cycle ({a}, lres, {dp}, relax, &p, nrelax, 0, depth());   
    foreach()
      if (fabs(res[]) > max)
	max = fabs(res[]);
      astats s = adapt_wavelet ({a}, (double []){cmax},minlevel=6,maxlevel=10,{a});
    Compute_B();
    printf("\n \n The Iteration Count : %d , The Max Residual : %g",i,max);
    printf("\n The number of refined grids : %d , The number of coarsened grids : %d, Total Grids : %d ",s.nf,s.nc,grid->tn);
  }
  
  foreach(){  
    if(x*x <=1 && y*y<=1)
      Error[] = a[] - (-pow(y,10)+5.*pow(y,8)-10.*pow(y,6)+10.*pow(y,4)-5.*pow(y,2)+1.)*(-pow(x,10)+5.*pow(x,8)-10.*pow(x,6)+10.*pow(x,4)-5.*pow(x,2)+1.);
    else
      Error[] = a[];
    }

  FILE *fp, *fp1, *fp2;
  fp = fopen("ErrorDomain.dat","w");
  fp1 = fopen("Residual.dat","w");
  fp2 = fopen("ErrorvsGrid.dat","a");       
  stats s2 = statsf(Error);
  foreach(){
    Error[] -= s2.sum/s2.volume;
    fprintf(fp,"%g %g %g \n",x,y,Error[]);
    fprintf(fp1,"%g %g %g \n",x,y,res[]);
    if (fabs(Error[]) > max) max = fabs(Error[]);
  }
  fprintf (fp2, "%g %g\n",log10(sqrt(grid->tn)),log10(normf(Error).max));
  fclose(fp);
  fclose(fp1);
  fclose(fp2);

}

int main (int argc, char ** argv)
{
  system("rm -f ErrorvsGrid.dat");
  for(double cmax=1E-03;cmax>=1E-06;cmax=cmax/10.)
    solve (cmax);
}

/**
## Results

~~~gnuplot Error Distribution for Residual operator
set output 'ErrorDistribution.png'
splot 'ErrorDomain.dat'
~~~

~~~gnuplot Error Distribution for Residual operator
set output 'ResidualDistribution.png'
splot 'Residual.dat'
~~~

~~~gnuplot Error Scaling of the Poisson 4th order solver on a heterogenous grid

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
