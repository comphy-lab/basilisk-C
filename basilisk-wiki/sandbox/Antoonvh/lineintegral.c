/**
#Evaluate a line integral on a 2D grid.
*/

#include "grid/quadtree.h"

#define func(a,b) (sin(a)*cos(b))
#define funcint(a,b) (sin(a)*(sin(Y0+L0)-sin(y-(Delta/2))))
scalar s[],ints[];
int startlevel = 2;

int main(){
  L0=2.*M_PI;
  X0=Y0=-L0/2;
  for (int i=1;i<7;i++){ 
    // Test to function on several grids
    init_grid(1<<(startlevel+i));
    refine(level<=startlevel+i && (x+y)<-1.);
    unrefine(level>=startlevel+i-1 && (x+y)>1.);
    
    // Define the function (2nd-order accurate)
    foreach()
      s[]=func(x,y);
    
    boundary({s}); //<- This is important
    // Integrate to the top
    foreach(){
      //Store height of this cell;
      double yt=y;
      double yh=y;
      double integral=0.;
      while(yh<(Y0+L0)){
        Point point = locate(x,yh);
        integral +=interpolate(s,x,y)*Delta;
        yh+=Delta;
      }
      Point point = locate(x,yt);
      ints[]=integral;
    }
    double err=0;
    foreach()
      err+=fabs(ints[]-funcint(x,y))*sq(Delta);
    static FILE * fp = NULL; if (!fp) fp = fopen("err.dat","w");
    fprintf(fp,"%d\t%g\n",i,err);
  }
}

/**
## Results
The function appears to be first-order accurate:

~~~gnuplot
set xr [ 0.5:6.5]
set logscale y
plot 'err.dat' u 1:2 ,\
 2**(-x)
~~~

*/
