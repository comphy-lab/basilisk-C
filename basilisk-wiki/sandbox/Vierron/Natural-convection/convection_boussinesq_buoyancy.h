#include "navier-stokes/centered.h"
#include "tracer.h"
#include "vof.h"
#include "diffusion.h"

/* Variables initialisation: T = temperature field, f = fraction field */
scalar T[], f[];
scalar * tracers = {T}, * interfaces = {f};
mgstats mgT;

double Ra, Pr, B;

event init (t=0) {
#if dimension == 2
  const face vector muc[] = {Pr/sqrt(Ra),Pr/sqrt(Ra)};
#endif
#if dimension == 3
  const face vector muc[] = {Pr/sqrt(Ra),Pr/sqrt(Ra),Pr/sqrt(Ra)};
#endif

  mu = muc;
  a = new face vector;
}

/* D = diffusion coefficient
   dt = time step
   T = Temperature field
   mgT = statistics of the Poisson solver */
event tracer_diffusion (i++) {
#if dimension == 2
  const face vector D[] = {1./sqrt(Ra), 1./sqrt(Ra)};
#endif
#if dimension == 3
  const face vector D[] = {1./sqrt(Ra), 1./sqrt(Ra), 1./sqrt(Ra)};
#endif
  mgT = diffusion (T, dt, D);
  boundary ({T});
}

/* Here, the acceleration term a of Navier-stokes/centered */
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] += Pr*((T[] + T[0,-1])/2.- B*(f[] + f[0,-1])/2.);
}
