/**
# The distance between two point within a cube. 
In order to test a double grid iterator implementation (i.e.nested) we aim to find the average distance between all points within a cube. The assumption here is that this number may be approximated using a finite number of points. The exact expression and a decimal expression is provided by [mathworld](http://mathworld.wolfram.com/HypercubeLinePicking.html), which will serve as a benchmark for our result. 
*/

#include "grid/octree.h"

int main(){
  FILE * fp = fopen("distance","w");
  int o = BGHOSTS+1;
  /**
  We check for an increasing number of points the distance to all points, using the cell centered coordinates. We increase the resolution upto 6 levels of refinement. Since we use a *rather naive* nested iterator, the scaling with increasing resolution is horrible (six dimensional).
  */
  for (int g=1;g<6;g++){
    int lev=g;
    init_grid(1<<lev);
    scalar xx[],yy[],zz[];
    int n=0;
    double dist=0;
    /**
    Since we are not allowed to nest two `foreach()` loops, it was decided to store the locations of the cell centeres in dedicated fields. I imagine persuiing a nested `foreach()` loop would be quite difficult, especially when working with a decomposed domain (e.g. when using MPI). 
    */
    foreach(){
      xx[]=x; yy[]=y; zz[]=z;
    }
      /**
   In order to determine the distances, we iterate over each grid cell and store the location as *xp,yp,zp*. 
    */
    foreach(){
      double xp=x; 
      double yp=y;
      double zp=z;
      /**
      For each cell, we loop over all cells using a coding friendly [Cartesian fashion](mortcscart.c) and determine the corresponding distance. Which is then added to the total distance variable *dist*. The iterator should also be able to iterate over tree-grid-cells with a variable resolution. 
      */
      for (int l =1;l<=depth();l++){
        int Nl= (1<<l);
        for (int i=o; i<Nl+o; i++){
          for (int j=o; j<Nl+o; j++){
            for (int k=o; k<Nl+o; k++){
              Point point;
              point.level=l; point.i=i; point.j=j; point.k=k;
              if(is_active(cell) && is_local(cell)){   //Does this cell exist on this pid?
                if(is_leaf(cell)){    //We only consider the leaf cells.
                  dist+=pow(sq(xp-xx[])+sq(yy[]-yp)+sq(zp-zz[]),0.5);
                  n++;
                }
              }
            }
          }
        }
      }
    }
    /**
    the result is written to a file, together with the error based on the analytical result. 
    */
    fprintf(fp,"%d\t%g %d\t%d\t%g\n",g,dist/(double)n,n,N,fabs(dist/(double)n-0.6617071822));
  }
}
/**
## Result

~~~gnuplot the Grid-based-approach results seems to converge towards the analytical solution with an increasing level of refinement
set xr [0.5:5.5]
set yr [0.5:0.7]
set ylabel 'Distance'
set xlabel 'Refinement Iteration'
plot   0.6617072 w lines lw 3 t 'Analytical results' ,\
       'distance' u 1:2 t 'Grid Approximated average distance'
~~~

We may focus a bit more on the convergence
~~~gnuplot The convergence rate with respect to number of evaluated pairs is not very impressive, bus as expected. 
set xr [32:5000000000] 
set yr [0.0001:0.2]
set logscale y
set logscale x
set ylabel 'Error'
set xlabel 'Number of pairs evaluated'
plot 'distance' u 3:5 t 'Error in the grid Approximated average distance' ,\
      x**(-0.333) t '2/6-power scaling'
~~~

Notice that for each refinement iteration the required effort for the algorithm increases with a factor of ($2^6$)=64, compared to the previous iteration. For the modestly-sized $32^3$ we have evaluated the distances over 1 billion pairs! also keeping in mind that the `pow(x,0.5)` function is not particularly cheap. If only we could come up with an adaptive-grid approach...
*/