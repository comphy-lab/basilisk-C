/**
# Direct Laplacian computation upto 4th Order - Multigrid

We compute the Laplacian using a 4th order formulation.
This exercise is performed to test the Error scalability of the formulation on a periodic domain for which
we use differently resolved grids & study how the maximum error behaves on the multigrid

*/

#include "grid/multigrid.h"
#include "Laplacian.h"
#include "utils.h"
#define BGHOSTS 2

void Solve(int depth){
  
   L0=1;
   init_grid(1<<depth);
   periodic(left);
   periodic(top);

   scalar A[], B[], Res[], dA[];
   foreach(){
       B[] = 0;
       A[] = cos(2.*pi*x)*cos(2.*pi*y);
    }
   boundary({A,B});
   Laplacian({A},{B});
   
   scalar Error[];
   foreach()
      Error[] = B[] + 2.*(sin(2.*pi*(x+Delta/2.))-sin(2.*pi*(x-Delta/2.)))*(sin(2.*pi*(y+Delta/2.))-sin(2.*pi*(y-Delta/2.)))/sq(Delta); 

   static FILE *fp;
   fp = fopen("ErrorvsGrid.dat","a");
   
   fprintf(fp,"%d %g %g %g %g \n",depth,pow(2,depth),normf(Error).max,log10(pow(2,depth)),log10(normf(Error).max));
   fclose(fp);
}

int main(){
  
 int depth;
 for(depth=3;depth<=10;depth++)
    Solve(depth);  
}

/**
## Results

~~~gnuplot Error Scaling of the Laplacian 4th Order on a Multigrid
set output 'ErrorScaling.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 4:5 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 4:5 w p t 'ErrorScaling', f(x) t title_f(a,b)
~~~

*/