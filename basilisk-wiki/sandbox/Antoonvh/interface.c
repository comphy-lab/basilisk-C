/**
## The area of an interface
Basilisk provides and extensive toolbox to do interface reconstruction. On this page we check the convergence properties for determining the surface area of an iso-surface countour of an analytical function.  
*/
#include "grid/octree.h"
#include "fractions.h"
#include "utils.h"
/**
More specifially, we try to numerically determine the surface area of a sphere with radius $R=0.5 $ a.u., such that the corresponding interface area should be $A=\pi\ \text{ a.u.}^2$.
We use a function that varies quadtratically in space and one that varies linearly with the radius.
*/
# define func (sq(x-xo) + sq(y-yo) + sq(z-zo) - sq(R)) 
# define func2 (pow(sq(x - xo) + sq(y - yo) + sq(z - zo), 0.5) - R)

double xo = 0.1, yo = M_PI/10.0;
double zo=-1./2.44, R = 0.5;
scalar f[], f2[];
int main(){
  L0 = 5.;
  X0 = Y0 = Z0 = -L0/2;
  init_grid(1 << 4);
  int i = 0;
  FILE * fp1 = fopen("interface.dat", "w");
  FILE * fp = popen("gfsview-batch3D interface.interface.gfv | ppm2gif --delay 200 > surface.gif","w");
  for (i=1; i<8; i++){
    double s = 1./pow(2., (double)i); 
    refine(fabs(func) < (s) && level < i + 4);
    fraction(f, func);
    fraction(f2, func2);
    int cells = 0;
    foreach()
      if (f[] > 0. && f[] < 1.)
	cells++;
    fprintf(fp1, "%d\t%d\t%g\t%g\n",i, cells, fabs(interface_area(f)-M_PI), fabs(interface_area(f2)-M_PI));
    output_gfs(fp);
    fprintf(fp, "Save stdout { format = PPM width = 600 height = 600}\n");
  }
  fclose(fp1);
}
/**
First we visually inspect the effect of the increasing spatial resolution on the corresoponding fidelity in the represenation of the spherical surface. 

![Visual inspectation of the spatial convergence. The surface appears to be better resolved with each refinement step.](interface/surface.gif)

Next we check if the number of interfactial cells that we have used does indeed scale as expected. 

~~~gnuplot 2D scaling for a 2D problem.
set xr [0.5:7.5]
set logscale y
set ylabel 'number of interfacial cells' 
set xlabel 'refinement iteration'
set key top left
plot   50*2**(2*x) with lines title '2D scaling' ,\
      "interface.dat" using 1:2 title 'used interfacial cells'
    
~~~

Now we check the error in the computed interfacial area: 

~~~gnuplot
set ylabel 'Error [a.u.^2]'
set key top right
plot "interface.dat" using 1:3 title 'Error in area f' ,\
     "interface.dat" using 1:4 title 'Error in area f2'
~~~

We see that it does not converge. This could be due to the fact that for the interface reconstruction from an analytical function, the *fractions()* function used (only) second-order-accurate interpolation to determine the 'edge' fraction of a cell's edge. Thus when we integrate over $O(n^2)$ cells, the second-order accuracy *per cell* results in $O(n^2) \times O(n^{-2}) = O(N^0)$, i.e zeroth order global accuracy!  
*/
