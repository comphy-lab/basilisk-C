/**
#Refine function original code - Higher order prolongation operator.
*/

#define BGHOSTS 2
#include "grid/bitree.h"
#include "utils.h"
#include "Adaptivity_New.h"

#if dimension == 1

double SolutionCompact (double x, double Delta){
  double xq1, xq2;
  xq1 = x - Delta/(2.*sqrt(3));
  xq2 = x + Delta/(2.*sqrt(3));
  if(sq(x) <=1)
    return ( 0.5 * ( pow((1-sq(xq1)),5) + pow((1-sq(xq2)),5) )  );
  else
    return (0);
}

double SolutionPeriodic (double x, double Delta){
 return ((2./(pi*Delta))*(cos(pi*x/2. - pi*Delta/4.)-cos(pi*x/2. + pi*Delta/4.)));
}

#endif



void convergence (int depth){

  L0 = 4;
  origin (-2);
  init_grid(1<<(depth-1));
  scalar s[];
  foreach()
     s[] = SolutionCompact (x,Delta);
  foreach_dimension(){
     s[left] = neumann(0);
     s[right] = neumann(0);
   }
  s.prolongation = refine_order5;
  boundary({s});

  refine((level < depth ) && (sq(x) <= 0.25) );								

  foreach_level(depth)
     s[] = Order5_Refine1D (point,s); 
  boundary({s});  

  scalar Error[];
  foreach()
     Error[] = s[] - SolutionCompact (x,Delta);

  FILE *fp = fopen("Error-Distribution.dat","w");
  foreach()
    fprintf(fp,"%g %g %g %g \n",x,s[],Error[],Delta);
  fclose(fp);

  fp = fopen("Error_Convergence.dat","a");
  fprintf(fp,"%g %g \n",log10(grid->tn),log10(fabs(normf(Error).max)) );
  fclose(fp);

  delete ({s,Error});
  free_grid();
}

int main(){
  ArrayInitializations();
  system ("rm -f Error_Convergence.dat");
  for (int depth = 7; depth <= 10 ; depth++)
     convergence (depth);
}

/**
## Results

~~~gnuplot Error distribution 1D - Full Domain
set key outside
set output '1DError-Order5.png'
set xlabel 'X'
set ylabel 'Error(X)'
set grid
plot 'Error-Distribution.dat' u 1:3 w p t 'Error-Distribution-Order5'
~~~

~~~gnuplot Solution construction 1D - Full Domain
set output '1DSolutionConstruction.png'
set key outside
set xlabel 'X'
set ylabel 'Delta(X)'
set grid
plot 'Error-Distribution.dat' u 1:2 w p t 'Solution-Distribution'
~~~

~~~gnuplot Delta-Distribution - Full Domain
set output '1D-Delta-Distribution.png'
set key outside
set xlabel 'X'
set ylabel 'Solution(X)'
set grid
plot 'Error-Distribution.dat' u 1:4 w p t 'Delta-Distribution'
~~~

~~~gnuplot Error Convergence of the prolongation operator
set key outside
set output 'Error-Convergence.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(L2Norm-Error)'
set grid
f(x) = a*x + b
fit f(x) 'Error_Convergence.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'Error_Convergence.dat' u 1:2 w p t 'Error-Convergence', f(x) t title_f(a,b)
~~~

*/