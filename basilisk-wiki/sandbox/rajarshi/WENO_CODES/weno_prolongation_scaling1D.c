/**
## Error Scaling of the Prolongation Operator !

In this code the error scaling of the WENO Prolongation operator has been demonstrated. 
*/

#include "grid/multigrid1D.h"
#include "utils.h"
#include "flux_construction.h"
#define BGHOSTS 2
#define dimension 1

double solution(double x, double Delta){
  return ( (2./(pi*Delta))*(cos(pi*x/2. - pi*Delta/4.) - cos(pi*x/2. + pi*Delta/4.)));
}

void Convergence(int depth){

   L0=4;
   origin(-2);
   init_grid(1<<depth);
   periodic(left);
   
   scalar s[],Error[];
   foreach_level_or_leaf(depth-1)
      s[] = solution(x,Delta);
   boundary_level({s},depth-1);
   
   s.prolongation = refine_weno;
   foreach_level(depth)
     s[] = WenoRefine1D(point,s); 
   boundary_level({s},depth);


   foreach()
      Error[] = s[] - solution(x,Delta);
   FILE *fp = fopen("Error-Scaling.dat","a");
   fprintf(fp,"%g %g \n",log10(pow(2,depth)),log10(normf(Error).max));
   fclose(fp);
  
  if(depth == 6){
   fp = fopen("Coarse.dat","w");
   foreach_level(depth-1)
     fprintf(fp,"%g %g \n",x,s[]);
   fclose(fp);

   fp = fopen("Fine.dat","w");
   foreach_level(depth)
     fprintf(fp,"%g %g \n",x,s[]);
   fclose(fp);
  }
}

int main(){
   system("rm -f Error-Scaling.dat");
   for(int depth=4; depth<=10; depth++)
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
*/