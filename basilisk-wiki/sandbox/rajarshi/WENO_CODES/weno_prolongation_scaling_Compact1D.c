/**
## Error Scaling of the Prolongation Operator

In this code the error scaling of the WENO Prolongation operator has been demonstrated. 
*/

#include "grid/multigrid1D.h"
#include "utils.h"
#include "flux_construction.h"
#define BGHOSTS 2

double solution(double xu, double xl, double Delta){
  return ((1./Delta)*((xu - (7./3.)*pow(xu,3) + (21./5.)*pow(xu,5) - 5.*pow(xu,7) + (35./9.)*pow(xu,9) - (21./11.)*pow(xu,11) + (7./13.)*pow(xu,13) - (1./15.)*pow(xu,15)) - (xl - (7./3.)*pow(xl,3) + (21./5.)*pow(xl,5) - 5.*pow(xl,7) + (35./9.)*pow(xl,9) - (21./11.)*pow(xl,11) + (7./13.)*pow(xl,13) - (1./15.)*pow(xl,15))));
}

void Convergence(int depth){

   L0=4;
   origin(-2);
   init_grid(1<<depth);
   periodic(left);
   
   scalar s[],Error[];
   foreach_level_or_leaf(depth-1){
     if(sq(x) <=1)
       s[] = solution(x+Delta/2.,x-Delta/2.,Delta);
     else
       s[] = 0.;
    }
   s[left] = neumann(0);
   s[right] = neumann(0);
   boundary_level({s},depth-1);
   
   s.prolongation = refine_weno;
   foreach_level(depth)
     s[] = WenoRefine(point,s); 
   boundary_level({s},depth);


   foreach(){
     if(sq(x)<=1)
        Error[] = s[] - solution(x+Delta/2.,x-Delta/2.,Delta);
     else
        Error[] = s[];
    }
   FILE *fp = fopen("Error-Scaling.dat","a");
   fprintf(fp,"%g %g \n",log10(pow(2,depth)),log10(normf(Error).max));
   fclose(fp);
   if(depth == 6){
     fp = fopen("Coarse.dat","w");
     foreach_level(depth-1)
        fprintf(fp,"%g %g \n",x,s[]);
     fclose(fp);

     fp = fopen("Fine.dat","w");
     foreach_level(depth){
       if(sq(x)<=1)
          fprintf(fp,"%g %g %g \n",x,s[],solution(x+Delta/2.,x-Delta/2.,Delta));
       else 
          fprintf(fp,"%g %g 0 \n",x,s[]);
       }
       fclose(fp);
    }

  fp = fopen("ErrorDistribution.dat","w");
  foreach()
    fprintf(fp,"%g %g \n",x,Error[]);
  fclose(fp);
}

int main(){
   system("rm -f Error-Scaling.dat");
   for(int depth=6; depth<=10; depth++)
      Convergence(depth);
}

/**
## Results

~~~gnuplot Error Scaling of the Prolongation Operator
set output 'ErrorScaling-Prolongationoperator.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'Error-Scaling.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'Error-Scaling.dat' u 1:2 w p t 'ErrorScaling', f(x) t title_f(a,b)
~~~

~~~gnuplot Comparing Coarse and Fine solutions
set output 'Coarse-Fine.png'
set xlabel 'x'
set ylabel 'f(x)'
set grid
plot 'Coarse.dat' u 1:2 w p t 'Coarse', 'Fine.dat' u 1:2 w p t 'Fine'
~~~

~~~gnuplot Comparing Analytical and Numerical solutions
set output 'AnalyticalVSNumerical.png'
set xlabel 'x'
set ylabel 'f(x)'
set grid
plot 'Fine.dat' u 1:2 w p t 'Numerical', 'Fine.dat' u 1:3 w p t 'Analytical'
~~~

~~~gnuplot Error distribution 1D - Full Domain
set output '1DError.png'
set xlabel 'X'
set ylabel 'Error(X) - Log Scale'
set grid
plot 'ErrorDistribution.dat' u 1:2 w p t 'Error-FullDomain'
~~~

*/
