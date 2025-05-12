/**
#General VOF-Riemann direction-splitted solver

We solve for the advection of a color function plus the conservative version of the
continuity, momentum and total energy equations
*/

#include "run.h"
#include "timestep.h"
#include "riemannsplitted.h"
#include "voftracerflux.h"

double dtmax;
scalar c[], * interfaces = {c}, * interfaces1 = {c};
scalar rho[], Etot[];
vector q[];
face vector uf[];

event defaults (i = 0) {
  foreach()
    rho[] = c[] = Etot[] = 1.;
  boundary ({rho,c, Etot});
}

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

event vof (i++) {

  face vector alpha[];
  foreach_face()
    alpha.x[] = 2./(rho[] + rho[-1]);
  boundary((scalar *){alpha});

  foreach_face()
    uf.x[] = alpha.x[]*(q.x[] + q.x[-1])/2.;
  boundary ((scalar *){uf});
  
  theta = 1.;
  foreach_dimension() {
    q.x.gradient = minmod2;
  }

  rho.gradient = minmod2;
  Etot.gradient = minmod2;

  c.tracers = {rho, Etot, q};

  vof_advection ({c}, i);

}