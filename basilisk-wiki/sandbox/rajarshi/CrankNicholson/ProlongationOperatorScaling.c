/**
## Error Scaling of the Prolongation Operator

In this code the error scaling of the prolongation operator has been demonstrated. There are currently
four different prolongation operators. The Bi-linear & Bi-quadratic ones are written in grids/multigrid.h
whereas the Bi-cubic and the Bi-quartic ones are written in Laplacian.h

The value of the parameter 'order' in the main function needs to be changed to observe the error scaling of the 
various schemes

*/

#include "grid/multigrid.h"
#include "utils.h"
#include "Laplacian.h"
#define BGHOSTS 2

double solution(double x, double y){
  return (cos(2.*pi*x)*cos(2.*pi*y));
}

void Convergence(int depth, int order){

   L0=1;
   origin(-0.5,-0.5);
   init_grid(1<<depth);
   foreach_dimension()
     periodic(left);
   scalar s[],Error[];
   foreach_level_or_leaf(depth-1)
      s[] = solution(x,y);
   boundary_level({s},depth-1);
   if (order==2){
     s.prolongation = refine_bilinear;
     foreach_level(depth)
      s[] = bilinear(point,s);
    }
   else if(order==3){
     s.prolongation = refine_biquadratic;
     foreach_level(depth)
      s[] = biquadratic(point,s);
    }
   else if (order==4){
     s.prolongation = refine_bicubic;
     foreach_level(depth)
      s[] = bicubic(point,s);
    }
   else if (order==5){
     s.prolongation = refine_biquartic;
     foreach_level(depth)
      s[] = biquartic(point,s);
     }
   else {
     printf("\n\n Enter a valid order (2, 3, 4 or 5)");
     abort();
     }
 
   boundary_level({s},depth);
   foreach()
      Error[] = s[] - solution(x,y);
   FILE *fp = fopen("Error-Scaling.dat","a");
   fprintf(fp,"%g %g \n",log10(pow(2,depth)),log10(normf(Error).max));
   fclose(fp);
}

int main(){
   remove("Error-Scaling.dat");
   int order = 5;
   for(int i=4;i<=10;i++)
      Convergence(i,order);
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

*/