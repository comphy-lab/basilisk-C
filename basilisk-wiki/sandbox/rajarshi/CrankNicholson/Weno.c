/**
# Weno
*/

#define BGHOSTS 2
#include "grid/multigrid1D.h"
#include "utils.h"

void Simple(scalar *al, face vector *bl){
 scalar a = al[0];
 face vector b = bl[0];
 foreach()
   a[] = 1.5;
 foreach_face()
   b.x[] = a[] - a[-1];
}

int main(){

 L0 = 4;
 origin(-2);
 init_grid(1<<5);
 scalar Z[];
 face vector X[];
 Simple({Z},{X});
 FILE *fp = fopen("Text.dat","w");
 foreach_face()
   fprintf(fp,"%g %g \n",x,X.x[]);
 fclose(fp);
}


/**
## Results

~~~gnuplot Fluctuations in Interpolation
set output 'Test.png'
set xlabel 'X'
set ylabel 'Z(X)'
set grid
plot 'Text.dat' u 1:2 w p t 'Test'
~~~

*/