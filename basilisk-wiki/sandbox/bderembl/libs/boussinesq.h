/**

   Boussinesq equations

*/

#include "tracer.h"
#include "diffusion.h"

face vector av[];
scalar b[];
scalar * tracers = {b};
(const) face vector kappa = zerof;

event init (t = 0) {
  a = av;
}

event acceleration (i++) {

#if dimension > 2
  coord grav = {0, 0, 1};
#else
  coord grav = {0, 1};
#endif

  foreach_face()
    av.x[] = grav.x*face_value(b, 0);
}


event tracer_diffusion (i++) {

  if (constant(kappa.x) != 0.) {
    mgstats mgb = diffusion (b, dt, kappa);
    boundary ({b});
  }
}
