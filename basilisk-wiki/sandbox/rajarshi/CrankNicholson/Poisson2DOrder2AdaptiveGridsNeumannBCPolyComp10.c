/**
# Poisson Solver 2nd Order Scheme on an Adaptive Grid

The second order Poisson scheme is used to solve a familiar compact function case with Neumann Boundary Conditions.
This case however uses adaptive grids which is done through the adapt_wavelet function.

*/

#include "grid/quadtree.h"
#include "poisson.h"
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
  b.prolongation = refine_bilinear;
  boundary ({b});

}

/**
This is the Solve function, which starts with setting the domain, and initializing an uniform grid with a given depth. Then a call is made to initialize the RHS.
Boundary conditions are applied on the function a (i.e. LHS, which needs to be computed) and a bilinear prolongation function is chosen from the header file Laplacian.h.
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
  a.prolongation = refine_bilinear;
  boundary({a,dp});


/**

The poisson problem is set up, and utilizing the functions residual and mg_solve in the header file poisson.h, a set of iterations are run to achieve convergent residuals. It is to be noted the outer loop carries out succesive stages of mesh adaptation ( the control criterion is set, such that further iterations stop once
the grids do not refine further. The inner loop controls the residual convergence on a given stage of adaptive refinement. Here the control criterion is set to a 
given tolerance of the residual. 

*/


  #define ITER 15
  #define NITER 20

  const scalar lambda[] = 0.;
  struct Poisson p;
  p.a = a; p.b = b; p.alpha = unityf; p.lambda = lambda;
  scalar * lres = {res};
  
  double max;
  int gridNo = 0;
  for (int i = 0; i < ITER && gridNo != grid->tn ; i++) {
    max=1;
    gridNo = grid->tn;
    for(int j = 0; j < NITER && max > 1E-09 ; j++){
       max=0;
       residual ({a}, {b}, lres, &p);
       mg_cycle ({a}, lres, {dp}, relax, &p, nrelax, 0, depth());
       foreach()
         if (fabs(res[]) > max)
	    max = fabs(res[]);
       printf("\n \n Adapatation Iteration : %d , The Iteration Count : %d , The Max Residual : %g",i,j,max);
       }
    astats s = adapt_wavelet ({a}, (double []){cmax},minlevel=depth,maxlevel=10,{a});
    Compute_B();
    printf("\n The number of refined grids : %d , The number of coarsened grids : %d, Total Grids : %ld ",s.nf,s.nc,grid->tn);
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
  fprintf (fp2, "%g %g\n",log10(grid->tn),log10(normf(Error).max));
  fclose(fp);
  fclose(fp1);
  fclose(fp2);

}

int main (int argc, char ** argv)
{
  system("rm -f ErrorvsGrid.dat");
  for(double cmax=1E-03; cmax>=1E-06; cmax=cmax/10.)
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