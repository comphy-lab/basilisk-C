/**
# Spark Model

For the ignition of fuel mixtures in combustion
simulations.
*/

#include "grid/multigrid.h"
#include "run.h"
#include "common.h"
#include "view.h"

/**
The spark must be applied on a temperature field.
*/

scalar T[];
scalar sgT[], rho2v[], cp2v[];
#include "spark.h"

int main (void) {

  /**
  We gather the spark properties in the following
  struct. We assign the temperature field, the
  position of the center of the spark, its diameter,
  the ignition time, the total durantion and the
  maximum temperature. 
  */

  spark.policy = 2; // SPARK_RAMP
  spark.T = T;
  spark.position = (coord){0.25, 0.25};
  spark.diameter = 4.*L0/(1 << 7);
  spark.time = 0.1;
  spark.duration = 0.1;
  spark.temperature = 3500.;

  /**
  We set a maximum delta t and we create the
  domain and the grid.
  */
  
  DT = 1.e-3;
  size (1.);
  origin (0., 0.);
  init_grid (1 << 7);
  run();
}

event init (i = 0) {
  /**
  The temperature is initially constant and at
  300K in the whole domain.
  */

  foreach()
    T[] = 300.;
  boundary({T});
}

event stability (i++,last) {
  dt = dtnext (DT);
}

/**
We write a video with the temperature field in time,
and we write the volume integral of the temperature
in the logfile.
![Temperature field](spark/movie.mp4)(width="800" height="600")
*/

event logfile (i++) {
  fprintf (stderr, "%f %f\n", t, statsf(T).sum);
}

event movie (t += 0.01; t <= 0.3) {
  clear ();
  view (tx = -0.5, ty = -0.5);
  squares ("T", min=300., max=spark.temperature, linear=true);
  save ("movie.mp4");
}

