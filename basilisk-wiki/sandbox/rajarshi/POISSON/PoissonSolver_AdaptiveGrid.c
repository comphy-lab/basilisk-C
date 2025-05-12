/**
# Poisson Solver on an Adaptive Grid
*/

#define ORDER 4

#include "grid/tree.h"

#if ORDER == 2
 #include "poisson.h"
 #define BGHOSTS 1
#elif ORDER == 4
 #include "PoissonO4.h"
 #define BGHOSTS 2
#endif

#include "utils.h"

#define dimension 2
#define Value 0
#define Choice 2   // Choice 1 : Periodic Functions. Choice 2 : Compact Functions

scalar a[], b[], res[], dp[], Error[];

/**
This is a function, which will be called whenever a re-computation of the RHS is required. Mostly this will be needed for initialization and after every adaptive refinement step.
*/

#if Choice == 1

double PoissonLHS (double x, double y, double Delta)
{
  return ( pow((1./(2.*pi*Delta)),dimension)*(cos(2.*pi*x - pi*Delta)-cos(2.*pi*x + pi*Delta))*(cos(2.*pi*y - pi*Delta)-cos(2.*pi*y + pi*Delta)) );
}

void PoissonRHS (){

 foreach() 
   b[] = -2.*pow((1./(Delta)),dimension)*(cos(2.*pi*x - pi*Delta)-cos(2.*pi*x + pi*Delta))*(cos(2.*pi*y - pi*Delta)-cos(2.*pi*y + pi*Delta));
 
 stats s = statsf(b);
 foreach()
    b[] -= s.sum/s.volume;
 boundary({b}); 
}

#elif Choice == 2

double PoissonLHS (double x, double y, double Delta){
  
  if(sq(x) <= 1 && sq(y) <= 1){
     double q = Delta/(2.*sqrt(3));
     return ( (pow((1-sq(x-q)),5)*pow((1-sq(y-q)),5) + pow((1-sq(x+q)),5)*pow((1-sq(y-q)),5) + pow((1-sq(x-q)),5)*pow((1-sq(y+q)),5) + pow((1-sq(x+q)),5)*pow((1-sq(y+q)),5) )/(4.) ); 
   }
  else
    return (0);
}


void PoissonRHS (){
   
  foreach() {
     if(sq(x) <= 1 && sq(y) <= 1){
         double q = Delta/(2.*sqrt(3));
         b[] = ( (  ( 80.*pow((1-sq(x-q)),3)*sq(x-q)*pow((1-sq(y-q)),5) - 10.*pow((1-sq(x-q)),4)*pow((1-sq(y-q)),5) + 80.*pow((1-sq(y-q)),3)*sq(y-q)*pow((1-sq(x-q)),5) - 10.*pow((1-sq(y-q)),4)*pow((1-sq(x-q)),5) ) + ( 80.*pow((1-sq(x+q)),3)*sq(x+q)*pow((1-sq(y-q)),5) - 10.*pow((1-sq(x+q)),4)*pow((1-sq(y-q)),5) + 80.*pow((1-sq(y-q)),3)*sq(y-q)*pow((1-sq(x+q)),5) - 10.*pow((1-sq(y-q)),4)*pow((1-sq(x+q)),5) ) + ( 80.*pow((1-sq(x-q)),3)*sq(x-q)*pow((1-sq(y+q)),5) - 10.*pow((1-sq(x-q)),4)*pow((1-sq(y+q)),5) + 80.*pow((1-sq(y+q)),3)*sq(y+q)*pow((1-sq(x-q)),5) - 10.*pow((1-sq(y+q)),4)*pow((1-sq(x-q)),5) ) + ( 80.*pow((1-sq(x+q)),3)*sq(x+q)*pow((1-sq(y+q)),5) - 10.*pow((1-sq(x+q)),4)*pow((1-sq(y+q)),5) + 80.*pow((1-sq(y+q)),3)*sq(y+q)*pow((1-sq(x+q)),5) - 10.*pow((1-sq(y+q)),4)*pow((1-sq(x+q)),5) ) )/(4.)   );
     }
    else
        b[] = 0;
  }
 
  stats s = statsf(b);
  foreach()
    b[] -= s.sum/s.volume;  
   
  boundary ({b});

}

#endif


/**
This is the Solve function, which starts with setting the domain, and initializing an uniform grid with a given depth. Then a call is made to initialize the RHS.
Boundary conditions are applied on the function a (i.e. LHS, which needs to be computed) and a biquartic prolongation function is chosen from the header file Laplacian.h.
The poisson problem is set up, and utilizing the functions residual and mg_solve in the header file poisson.h, a set of iterations are run to achieve a convergent residual.
These calculations are all on the uniform Multigrid.
*/

void solve (double cmax)
{ 
  
 int nrelax = 4;
 clock_t start, end;
 start = clock(); 

 #if Choice == 1

  origin (-1,-1);
  L0 = 2;
  init_grid(1<<6);
  foreach_dimension()
    periodic(left);

 #elif Choice == 2

  origin (-2.5,-2.5);
  L0 = 5;
  init_grid(1<<6); 
  foreach_dimension(){
    a[left] = neumann(Value);
    a[right] = neumann(Value);
    b[left] = neumann(Value);
    b[right] = neumann(Value);
  }

 #endif

 #if ORDER == 4
  for(scalar x in {a,b}){
     x.refine = refine_order5;
     x.prolongation = refine_order5;
  }
 #endif

  PoissonRHS();
  foreach()
    a[] = 0.;
  boundary({a});


  const scalar lambda[] = 0.;
  struct Poisson p;
  p.a = a; p.b = b; p.alpha = unityf; p.lambda = lambda;
  scalar * lres = {res};

  double maxnow, maxprev, tolerance;
  astats s;
  s.nf = s.nc = 2;

  while(s.nf >=1 || s.nc >=1){
    maxnow=1;
    maxprev=0;
    tolerance = 1;    
    while(tolerance >= 1E-14){
       residual ({a}, {b}, lres, &p);
       mg_cycle ({a}, lres, {dp}, relax, &p, nrelax, 2, depth());       
       maxnow=0;
       foreach()
         if (fabs(res[]) > maxnow)
	    maxnow = fabs(res[]);
       tolerance = fabs(maxnow-maxprev);
       maxprev = maxnow;
    }
    
    s = adapt_wavelet ({a}, (double []){cmax},minlevel=5,maxlevel=18,{a});
    PoissonRHS();
   }
  
  foreach()
    Error[] = a[] - PoissonLHS(x,y,Delta);
  
  FILE *fp2 = fopen("ErrorvsGrid.dat","a");       
  stats s2 = statsf(Error);
  foreach()
    Error[] -= s2.sum/s2.volume;
  end = clock();
  
  if(cmax == 4E-05){
    FILE *fp, *fp1;
    fp = fopen("ErrorDomain.dat","w");
    fp1 = fopen("Residual.dat","w");  
    foreach(){
      fprintf(fp,"%g %g %g %g %g \n",x,y,a[],b[],Error[]);
      fprintf(fp1,"%g %g %g \n",x,y,res[]);
    }
    fclose(fp);
    fclose(fp1);    
  }

  fprintf (fp2, "%g %g\n",log10(sqrt(grid->tn)),log10(normf(Error).max));
  fclose(fp2);
  
  FILE *fp3 = fopen("TimeComputing.dat","a");
  fprintf (fp3, "%g %g\n",(end-start)/(double)(CLOCKS_PER_SEC),normf(Error).max );
  fclose(fp3);


}

int main (int argc, char ** argv)
{
  system("rm -f ErrorvsGrid.dat");
  system("rm -f TimeComputing.dat");
  #if ORDER == 4
    Prolongation_Weight_Initialization();
  #endif
  for(double cmax=1E-03 ; cmax>=1E-07 ; cmax=cmax/5.) 
    solve (cmax);
}

/**
## Results

~~~gnuplot Error Distribution for Residual operator
set output 'ErrorDistribution.png'
splot 'ErrorDomain.dat' u 1:2:5 w p t 'Error'
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

~~~gnuplot Computation time required

set output 'TimeComputation.png'
set logscale xy
set xlabel 'Time'
set ylabel 'Error'
set grid
plot 'TimeComputing.dat' u 1:2 w lp t 'Time(Poisson-O4)'
~~~
*/