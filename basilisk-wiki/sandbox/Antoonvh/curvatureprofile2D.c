/**
# Curvature profile of a 1D interface
If we have a nicely resolved interface, we could aim to to find it's curvature along a the interface. 

## set-up
Include the required header files. 
*/

#include "grid/quadtree.h"
#include "fractions.h"
#include "curvature.h"
#include "utils.h"
#include "interface_iterator.h"

/**
Define functions that will help create example interfaces. 

Shape 1 has a functional relation between y and x:

$$ y = \frac{1}{x+L0/2+0.01}+\mathrm{sin}(3x)$$

Shape 2 is a non-centred oval:

$$ (1.1(x-0.01))^2 + (0.9y)^2 = 1 $$

*/

#define func (1./(x+L0/2+0.01)+sin(3*x)+(-y)) // y=f(x) shape
#define func2 (sq(1.1*(x-0.01))+sq(0.9*y)-1.) // Bubble-like shape

int main(){
  /** 
  Set-up the domain parameters and the all-important starting location. The interfacial cell that is not at the domain boundary, whos centred location is *closest* to this *starting location* will be the starting cell. 
  
  The initial 'take-off' direction of the iterator will be away from this starting coordinate, i.e. initialized in *xyn[2]*. We choose it to be the top-left corner, that is not connected to an interfacial cell itself! 
  
  The algorithm tries to continue in this direction, following the interface. However, 'bounce backs' and other unexpected behaviour can occur for poorly resolved interfaces. 
  */
  L0=2*M_PI;
  X0=Y0=-L0/2.;
  double xyn[2];
  xyn[0] = X0;
  xyn[1] = Y0+L0;
  init_grid(64);
  scalar f[],tag[];
  f.refine=f.prolongation=fraction_refine;
  
  /**
  ## First shape
  */
  fraction(f,func);
  boundary({f});
  while(adapt_wavelet({f},(double[]){0.005},10).nf > 100){ 
    fraction(f,func);
    boundary({f});
  }
  curvature(f,tag);
  FILE * fp = fopen("result1","w");
  loop_interfacial_cells(fp,f,tag,xyn);
  
  /**
  ## Second shape
  We challenge the function by not adapting the grid with respect to this interface shape. 
  */
  fraction(f,func2);
  boundary({f});
  curvature(f,tag);
  fp = fopen("result2","w");
  loop_interfacial_cells(fp,f,tag,xyn);
}

/**
## Results
First we check if the interface is iterated in a sensible sequence by connecting the subsequent interface-centre locations with lines. 

~~~gnuplot The interface, reconstructed with a connected line. 
set xr [-3.14:3.14]
set yr [-3.14:3.14]
set xlabel('x')
set ylabel('y')
set size ratio -1
set key box bottom right
plot "result1" using 2:3 with lines title 'The functional interface' lw 3,\
     "result2" using 2:3 with lines title 'The Bubble-like interface' lw 3
~~~

That looks ok,complete(!) and no wierd jumps etc. Next we check if we can spot the minima and maxima in the curvature of the red curve, visually:

~~~gnuplot Curvature along the interface
set xr [0:16]
set yr [-20:12]
set xlabel('Interfacial distance from start')
set ylabel('Curvature [a.u.], x and y coordinates')
set size square

plot "result1" using 4:5 with lines lw 3 title 'Curvature' ,\
"result1" using 4:2 with lines title 'x coordinate' ,\
"result1" using 4:3 with lines title 'y coordinate'
~~~

That looks OK-ish. It would be relatively straightforward to find local extreme values and associate them with the centred coordinates of the individual reconstructed interfaces.

For the bubble-like interface we check the correlation between x-coordinate and curvature:

~~~gnuplot Curvature vs x-coordinate
set xr [-1:1]
set yr [-2:-0.5]
set xlabel('x')
set ylabel('Curvature [a.u.]')
set size square
plot "result2" using 2:5 
~~~

very usefull. 

## Note
The algorithm has a lot of quircks, please read the code under "interface_iterator.h" carfully, so you can understand when you obtain rubbish. The function is MPI compatible, but it is not parallelized! Actually, there is only additional overhead (all-to-all communications) when using multiple threads.   

I have only tested a few interfaces and errors typically occured. On one hand, the algorithm is not very robust, on the otherhand, for poorly resolved interfaces; shit in equals shit out. 

*/