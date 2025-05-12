/**
# Resolution of the 2D diffusion equation 
Implicit resolution of the diffusion equation

$$\frac{\partial T}{\partial t}= \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2}$$ 

using the Poison solver (see diffusion.h)

*/


#include "diffusion.h"
#include "run.h"

#define EPS 0.1

scalar f[];

const face vector D[] = { 1. , 1. };

/** Analytical solution 

$$ f = e^{r^2/4. t}/(4. \pi t ) $$

*/

	double solution (double x, double y, double t)
	{
  	return  exp( -1. * (sq(x) + sq(y))/ (4. * t))/ (4. * 	pi * t) ;
	}

/** We use this function to initialize the computation (we use t=0.1).
*/

	event init (t = 0)
	{
  	foreach()
    	f[] = solution(x,y,0.1);
  	boundary ({f});
}

/** Running */
 
	event running ( i++ )
	{
          dt = dtnext (t, 0.01);
          diffusion (f,dt,D);
          boundary ({f});
        }

/** Output
every 0.1 
*/

	event print ( t = 0.1 ; t += 0.1 ; t <= 1. )
	{
  	double shift = 0.1 ;
  	// For y=0
  	for (double x = -L0/2 ; x <= L0/2; x += L0/200.)
    		{
      	printf ("%f %f %f \n", x, interpolate (f, x, 	0.0),solution(x,0.0,t+shift));
    	}
  	printf ("\n\n");
	}

/** Main program
Compile : qcc -lm diffusion.c -o diffusion

Run : diffusion > difussion.dat

Plotting : plot "difussion.dat" i 4 u 1:2 t "Computation" w p,"" i 4 u 1:3 t "Theory" w l
*/
	int main() {
  	// Lenght
  	L0 = 10.;
  	// coordinates of lower-left corner
  	X0 = Y0 = -L0/2;
  	//
  	N = 128*2 ;

  	run();
	}

/** 
![Comparison simulation-theory](/sandbox/diffusion.png)
        

*/