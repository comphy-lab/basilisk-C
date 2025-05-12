#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar T[];
scalar * tracers = {T};
mgstats mgT;

double Ra, Pr, Re;

event properties (i++) {

        mu = new face vector;
  	face vector muc = mu;
	foreach_face()
		muc.x[] = fm.x[]/Re;
		
	boundary ((scalar*){muc});
	mu = muc;
}


event acceleration (i++) {
	a = new face vector;
	face vector av = a;
	foreach_face(y)
		av.y[] += (Ra*Pr/sq(Re))*(T[] + T[0,-1])/2.;
	a = av;
}

event tracer_diffusion (i=0) {

	face vector D = new face vector;
	foreach_face()
		D.x[] = fm.x[]/(Re*Pr);
	boundary ((scalar*){D});
	
	mgT = diffusion (T, dt, D);
  	boundary ({T});
}