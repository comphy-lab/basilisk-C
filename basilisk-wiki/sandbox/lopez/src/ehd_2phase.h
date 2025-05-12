/**
# EHD of two-phase systems

This files serves to setup more easily EHD problems with interfaces.

The interface is tracked with a Volume-Of-Fluid method we assume that
the fluid 1 is of permittivity *epsilon1* and conductivity *cond1* and
the fluid 2 is a fluid of permittivity *epsilon2* and conductivity *cond2*. */

#include "two-phase.h"
#include "implicit.h"
#include "stress.h"
#ifdef FRACFACE
#include "fracface.h"
#endif

double  epsilon1 = 1, epsilon2 = 1, cond1 = 0, cond2 = 0.;

#ifndef perm
# define perm(f) (clamp(f,0.,1.)*(epsilon1 - epsilon2) + epsilon2)
#endif
#ifndef cond
# define cond(f) (clamp(f,0.,1.)*(cond1 - cond2) + cond2)
#endif

face vector epsilonv[];
scalar * tracers = {rhoe};

event defaults (i=0) {
  epsilon = epsilonv;
  if (cond1 || cond2)
    K = new face vector;

  f.tracers = f.tracers ? list_concat (f.tracers, tracers) : tracers;
}

event properties (i++) {
#ifdef FRACFACE
  face vector cf[];
  face_fraction (f, cf);
  foreach_face(){
    double ff = (f[] + f[-1])/2.;
    epsilon.x[] = perm(ff)*fm.x[];
    if (cond1 || cond2)
      K.x[] = cond(cf.x[])*fm.x[];
  }
#else
  foreach_face(){
    double ff = (f[] + f[-1])/2.;
    epsilon.x[] = perm(ff)*fm.x[];
    if (cond1 || cond2)
      K.x[] = cond(ff)*fm.x[];
  }
#endif  
}
