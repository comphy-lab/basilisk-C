/**
# Source of a river

The case studied here is a topography descending towards the y decreasing with a
parabolic bed of which one digs a little more the center. We use the function 
eta_b to impose a flow on the top edge. We study the water elevation to see if 
it is regular. We check that the imposed inflow is equal to the one entering the 
domain. 
*/

#include "saint-venant.h"
#include "discharge.h"

#define LEVEL 8

int main()
{
  size (10.);
  origin (- L0/2., - L0/2.);
  G = 9.81;
  N = 1 << LEVEL;
  run();
}



/**
The flow rate varies in time and is set by computing the the elevation
$\eta_s$ of the water surface necessary to match this flow rate. */

double etas;
double inflow_imp;

event inflow (i++) { 
	if (t < 1)
		inflow_imp = 0.5*t;
	else
		inflow_imp = 0.5;
	etas = eta_b (inflow_imp, top);
	h[top] = max (etas - zb[], 0.);
	eta[top] = max (etas - zb[], 0.) + zb[];
}

/**
## Initial conditions


*/

event init (i = 0)
{  
  
  DT = 1e-2;

  foreach(){
  	  if (x*x < 1)
  	  	  zb[] = 5*(x+1)*(x-1)+(x*x*0.01 + y)/2.;
 		  	  	 
  	  else
  	  	  zb[] = (x*x*0.01 + y)/2.;
  }
  boundary ({zb});
  
  output_ppm (zb,file = "topo.ppm");
  ![Evolution of the water level](riverbug/topo.ppm)
}

/**
## Outputs

We compute the time-derivative of the total water volume (i.e. the net
flow rate), and make a GIF movie. */

event logfile (i++; t <= 2.) {
  static double volo = 0., to = 0.;
  double vol = statsf(h).sum;
  if (i > 0)
    fprintf (ferr, "%g %.6f %.6f %g %g\n", t, (vol - volo)/(t - to), inflow_imp,vol, etas);
  volo = vol, to = t;
}

event output (i += 5) {
  output_ppm (h, min = 0, max = 0.1, file = "riverbug.gif");
}


/**
## Results
The inflow is different from the one imposed.

~~~gnuplot Evolution of the flow rate
set key top left
set xlabel 'Time'
set ylabel 'Flow rate'
plot './log' u 1:2 w l t 'obtained', \
	 './log' u 1:3 w l t 'imposed'
~~~

The flow is irregular 

![Evolution of the water level](riverbug/riverbug.gif)

## See also


*/