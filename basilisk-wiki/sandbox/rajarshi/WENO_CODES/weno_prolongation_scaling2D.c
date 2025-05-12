/**
## Error Scaling of the Prolongation Operator

In this code the error scaling of the WENO Prolongation operator has been demonstrated for two dimensions. 

*/

#include "grid/multigrid.h"
#include "utils.h"
#include "flux_construction.h"
#define dimension 2
#define BGHOSTS 2

double solution(double x, double y, double Delta){
  return ( (4./(pi*pi*Delta*Delta))*(cos(pi*x/2. - pi*Delta/4.) - cos(pi*x/2. + pi*Delta/4.)) * ( cos(pi*y/2. - pi*Delta/4.) - cos(pi*y/2. + pi*Delta/4.)) );
}

void Convergence(int depth){

   L0=4;
   origin(-2,-2);
   init_grid(1<<depth);
   periodic(left);
   periodic(top);

   scalar s[],Error[];
   foreach_level_or_leaf(depth-1)
      s[] = solution(x,y,Delta);
   boundary_level({s},depth-1);
   
   s.prolongation = refine_weno;
   foreach_level(depth)
     s[] = WenoRefine2D(point,s); 
   boundary_level({s},depth);


   foreach()
      Error[] = s[] - solution(x,y,Delta);

   FILE *fp = fopen("Error-Scaling.dat","a");
   fprintf(fp,"%g %g \n",log10(pow(2,depth)),log10(normf(Error).max));
   fclose(fp);

   if (depth == 5){

   fp = fopen("Coarse.dat","w");
   foreach_level(depth-1)
     fprintf(fp,"%g %g %g \n",x,y,s[]);
   fclose(fp);

   fp = fopen("Fine.dat","w");
   foreach_level(depth)
     fprintf(fp,"%g %g %g \n",x,y,s[]);
   fclose(fp);
   }

   fp = fopen("Error-Distribution.dat","w");
   foreach()
     fprintf(fp,"%g %g %g \n",x,y,Error[]);
   fclose(fp);
  
}


int main(){
   system("rm -f Error-Scaling.dat");
   for(int depth=3; depth<=10; depth++)
      Convergence(depth);
}

/**
## Results

~~~gnuplot Comparing Coarse and Fine solutions
set output 'Coarse-Fine.png'
set xlabel 'x'
set ylabel 'y'
set zlabel 'f(x,y)'
splot 'Coarse.dat' u 1:2:3 w p t 'Coarse-Solution', 'Fine.dat' u 1:2:3 w p t 'Fine-Solution'
~~~

~~~gnuplot Observing Error Distribution
set output 'Error-Distribution.png'
set xlabel 'x'
set ylabel 'y'
set zlabel 'E(x,y)'
splot 'Error-Distribution.dat' u 1:2:3 w p t 'Error-Distribution'
~~~

~~~gnuplot Error Convergence of the Prolongation Operator
set output 'ErrorConvergence-Prolongationoperator.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'Error-Scaling.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'Error-Scaling.dat' u 1:2 w p t 'Error-Convergence', f(x) t title_f(a,b)
~~~
*/

