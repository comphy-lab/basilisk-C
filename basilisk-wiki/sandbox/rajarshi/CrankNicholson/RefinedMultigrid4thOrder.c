/**
## POISSON SOLVER ON A SUCCESSIVELY REFINED GRID

The domain is made up of circularly refined grids starting from the centre. The Central core has the highest grid density, followed by
a level of lesser grid densities as we reach the edge of the domains.
In this case 3 layers of refinement are used for each domain.
The error scaling is studies as the grids for different levels (ranging from 8 to 10) are used for the central circle 
*/

#include "utils.h"
#include "Laplacian.h"
#include "grid/tree.h"
#include "poisson-helmholtz4th.h"
#define BGHOSTS 2
#define value 0

scalar a[], b[], res[], dp[];

double solution (double x, double y)
{
  return cos(2.*pi*x)*cos(2.*pi*y);
}

/** 
 The solve function starts with a given grid density at the centre (signified by the argument depth) and there are two lower levels 
 of refinement (depth-1 & depth-2) for the outer domains.
*/

void solve (int depth)
{ 
  origin (-0.5, -0.5, -0.5);
  int nrelax = 4;
  init_grid(1);
  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));
  
  foreach() {
    b[] = -8*pi*pi*cos(2.*pi*x)*cos(2.*pi*y);
    a[] = 0;
  }

  stats s = statsf(b);
  foreach()
    b[] -= s.sum/s.volume;  
 
  foreach_dimension()
     periodic(left);

  a.prolongation = refine_biquartic;
  boundary ({a,b});

  #define NITER 30
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  const scalar lambda[] = 0.;
  struct Poisson p;
  p.a = a; p.b = b; p.alpha = unityf; p.lambda = lambda;
  scalar * lres = {res};
  
  for (int i = 0; i < NITER; i++) {
    residual ({a}, {b}, lres, &p);
    mg_cycle ({a}, lres, {dp}, relax, &p, nrelax, 0, depth());
    double max = 0.;
    foreach(reduction(max:max))
      if (fabs(res[]) > max)
	max = fabs(res[]);
    iter[i] = clock();
    maxres[i] = max;
  }

/** 
   The Errors on the domain are calculated and written to files. While the Max Error for different levels of central refinement are written to another file and they are then compared to look at the error scaling of the finite volume solver.

*/
         
  double max = 0;
  static FILE *fp, *fp1, *fp2;
  fp = fopen("ErrorDomain.dat","w");
  fp1 = fopen("Residual.dat","w");
  fp2 = fopen("ErrorvsGrid.dat","a");       
  foreach(reduction(max:max)) {
    double e = a[] - solution(x, y);
    fprintf(fp,"%g %g %g \n",x,y,e);
    fprintf(fp1,"%g %g %g \n",x,y,res[]);
    if (fabs(e) > max) max = fabs(e);
    // printf ("%g %g %g %g %g %g\n", x, y, a[], b[], res[], e);
  }
  fprintf (fp2, "%g %g %g %g\n", pow(2,depth), max, log10(pow(2,depth)),log10(max));
  fclose(fp);
  fclose(fp1);
  fclose(fp2);

}

int main (int argc, char ** argv)
{
  system("rm -f ErrorvsGrid.dat");
  for (int depth = 6; depth <= 9; depth++)
    solve (depth);
}

/**
## Results

~~~gnuplot Error Distribution for Residual operator
set output 'ErrorDistribution.png
splot 'ErrorDomain.dat'
~~~

~~~gnuplot Error Distribution for Residual operator
set output 'ResidualDistribution.png
splot 'Residual.dat'
~~~

~~~gnuplot Error Scaling of the Poisson 4th order solver on a heterogenous grid

set output 'ErrorScaling.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 3:4 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 3:4 w p t 'ErrorScaling', f(x) t title_f(a,b)
~~~

*/

