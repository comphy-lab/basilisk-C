#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar f[]; // This is the density difference
scalar * tracers = {f};
mgstats mgT;

double Gr, Sc ;

event init (t=0) {
#if dimension == 2
  const face vector muc[] = {1./sqrt(Gr),1./sqrt(Gr)};
#endif

  mu = muc;
  a = new face vector;
}

event tracer_diffusion (i++) {
#if dimension == 2
	const face vector D[] = {1./sqrt(Gr*Sc*Sc), 1./sqrt(Gr*Sc*Sc)};
#endif

  mgT = diffusion (f, dt, D);
  boundary ({f});
}


// Overload the acceleration event

event acceleration (i++) {
	face vector av = a;
	foreach_face(y)
		av.y[] -= (f[] + f[0,-1])/2.;;
}