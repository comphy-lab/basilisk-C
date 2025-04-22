/**
# Source of a river

In this example we impose a variable flow rate for a "source" located
inside the computation domain. */

#include "saint-venant.h"
#include "discharge.h"

/**
The domain is 10 metres squared, centered on the origin. Time is in
seconds. */

#define LEVEL 8

int main()
{
  size (10.);
  origin (- L0/2., - L0/2.);
  G = 9.81;
  N = 1 << LEVEL;
  
  /**
  We start with a dry riverbed, so that the problem does not have a
  natural timescale the Saint-Venant solver can use. We set a maximum
  timestep to set this timescale. */
  
  DT = 1e-2;
  
  run();
}

/**
## Boundary conditions

We create a new boundary for the source, with a Neumann condition for
the normal velocity (i.e. an inflow). */

bid source;
u.n[source] = neumann(0);

/**
The flow rate varies in time and is set by computing the the elevation
$\eta_s$ of the water surface necessary to match this flow rate. */

double etas;

event inflow (i++) {
  etas = eta_b (0.1*(1.1 - cos(4.*pi*t)), source);
  h[source] = max (etas - zb[], 0.);
  eta[source] = max (etas - zb[], 0.) + zb[];
}

/**
## Initial conditions

The river bed is a single valley. The source is created by
[masking](/Basilisk C#complex-domains) and is a narrow slot located at
the head of the valley. */

event init (i = 0)
{
  mask (fabs(x) < 0.5 && fabs(y - 3.5) < Delta/2. ? source : none);
  foreach()
    zb[] = (- cos(x) + y)/2.;
}

/**
## Outputs

We compute the time-derivative of the total water volume (i.e. the net
flow rate), and make a GIF movie. */

event logfile (i++; t <= 2.) {
  static double volo = 0., to = 0.;
  double vol = statsf(h).sum;
  if (i > 0)
    fprintf (stderr, "%g %.6f %g %g\n", t, (vol - volo)/(t - to), vol, etas);
  volo = vol, to = t;
}

event output (i += 5) {
  output_ppm (h, min = 0, max = 0.05, file = "source.gif");
}

/**
## Results

~~~gnuplot Evolution of the flow rate
set key top left
set xlabel 'Time'
set ylabel 'Flow rate'
plot './log' u 1:2 w l t 'obtained', 0.1*(1.1-cos(4.*pi*x)) w l t 'imposed'
~~~

![Evolution of the water level](source/source.gif)

## See also

* [Flow rates for multiple rivers](multiriverinflow.c)
*/
