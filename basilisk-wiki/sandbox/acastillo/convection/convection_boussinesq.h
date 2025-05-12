#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar T[];
scalar * tracers = {T};
mgstats mgT;

double Ra, Pr;

event init (t=0) {

	a = new face vector;
  mu = new face vector;

  face vector muc = mu;
  foreach_face()
    muc.x[] = fm.x[]*Pr/sqrt(Ra);
  boundary ((scalar*){muc});
}

event tracer_diffusion (i++) {

  face vector D[];
  foreach_face()
    D.x[] = fm.x[]*1./sqrt(Ra);
  boundary ((scalar*){D});
#if EMBED
	mgT = diffusion (T, dt, D, theta=cm);
#else
	mgT = diffusion (T, dt, D);
#endif
  boundary ({T});
}

event acceleration (i++) {
	face vector av = a;
	foreach_face(y)
		av.y[] += Pr*(T[] + T[0,-1])/2.;
}
