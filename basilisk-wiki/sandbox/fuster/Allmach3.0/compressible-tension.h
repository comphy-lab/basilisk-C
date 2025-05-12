/**
#Surface tension effects using the compressible solver

*/

#include "tension.h"

event end_timestep (i++)
{
  face vector pf[];

  scalar skappa[];
  foreach ()
    skappa[] = 0.;
  if (f.sigma > 0.)
    curvature (f, skappa, f.sigma, add = true);
  boundary({skappa});

  foreach_face () {
    pf.x[] = 0.;
    if ((f[1] + f[])/2. > 0. && (f[1] + f[])/2. < 1. ) {
      if (f[1] > f[])
        pf.x[] = (1. - f[1])*skappa[1];
      else
        pf.x[] = (1. - f[])*skappa[];
    }
  }
  boundary((scalar *){pf});

  foreach () {
    double div1 = 0.;
    foreach_dimension () 
      div1 += (pf.x[1]*uf.x[1] - pf.x[]*uf.x[]);

    fE1[] -= f[]*div1/Delta*dt/cm[];
  }

  foreach_face () {
    pf.x[] = 0.;
    if ((f[1] + f[])/2. > 0. && (f[1] + f[])/2. < 1. ) {
      if (f[1] > f[])
        pf.x[] = - f[]*skappa[];
      else
        pf.x[] = - f[1] *skappa[1];
    }
  }
  boundary((scalar *){pf});

  foreach () {
    double div1 = 0.;
    foreach_dimension () 
      div1 += (pf.x[1]*uf.x[1] - pf.x[]*uf.x[]);

    fE2[] -= (1. - f[])*div1/Delta*dt/cm[];
  }

  boundary({fE1,fE2});
}
