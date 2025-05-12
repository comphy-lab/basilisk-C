#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar T[];
scalar * tracers = {T};
mgstats mgT;

double Ra, Pr;

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

event acceleration (i++) {
	face vector av = a;
	foreach_face(y)
		av.y[] += Pr*(T[] + T[0,-1])/2.;
}