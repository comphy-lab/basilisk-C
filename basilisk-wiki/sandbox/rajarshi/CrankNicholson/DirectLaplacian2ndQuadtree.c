/**
## Compuation of Direct Laplacian Upto 2nd Order on a non homogenous grid (Quadtree)

The domain is made up of refined grid patches starting from the centre. The Central core has the highest grid density, followed by a level of lesser grid densities as we reach the edge of the domains. In this case 2 layers of refinement are used for each domain. The error scaling is studied as the grids for different levels (ranging from 6 to 10) are used for the central circle.
*/

#include "Laplacian.h"
#include "utils.h"
#define BGHOSTS 2

void Solve(int depth){
  
   L0=1;
   origin(-0.5,-0.5);
   init_grid(1);
   refine (level < depth - 1 || (sqrt(x*x)<0.125 && level < depth) );
   periodic(left);
   periodic(top);

   scalar A[], B[], Res[], dA[];
   foreach(){
       B[] = 0;
       A[] = cos(2.*pi*x);
    }
   A.prolongation = refine_bilinear;
   boundary({A,B});
   Laplacian2nd({A},{B});
   scalar Error[];
   foreach()
      Error[] = B[] + 4.*pi*pi*cos(2.*pi*x);
   FILE *fp,*fp2;
   fp = fopen("ErrorvsGrid.dat","w");
   fprintf(fp,"%d %g %g %g %g \n",depth,pow(2,depth),normf(Error).max,log10(pow(2,depth)),log10(normf(Error).max));
   fclose(fp);
   fp2 = fopen("ErrorDistribution.dat","w");
   foreach()
     fprintf(fp2,"%g %g %g \n",x,y,Error[]);
   fclose(fp2);
}

int main(){
 remove("ErrorvsGrid.dat");
 int depth;
 for(depth=7;depth<=7;depth++)
    Solve(depth);  
}

/**
## Results

~~~gnuplot Error Distribution of Laplacian
set output 'ErrorDistribution.png'
splot 'ErrorDistribution.dat'
~~~
*/