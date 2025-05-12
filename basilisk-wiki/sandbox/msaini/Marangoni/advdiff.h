#include "bcg.h"
#include "diffusion.h"

extern face vector uf;
extern double dt;
extern double rho1,rho2;
extern scalar f;

double Diff_C1 = 0, Diff_C2 = 1.;
scalar PhiC[];

#ifndef Diffvv
# define Diffvv(f) (clamp(f,0.,1.)*(Diff_C1 - Diff_C2) + Diff_C2)
#endif


#if TREE
event defaults (i = 0)
{
    PhiC.refine  = refine_bilinear;
    PhiC.restriction = restriction_volume_average;
    //    PhiC.gradient = p.gradient;
    PhiC.gradient = zero;
    PhiC.dirty = true;
}
#endif // TREE

event tracer_advection(i++){
  advection ((scalar *){PhiC}, uf, dt);
}

face vector Diffv[];

event tracer_diffusion(i++){
  foreach_face(){
    double ff = (f[] + f[-1])/2.;
    Diffv.x[] = fm.x[]*Diffvv(ff);
  }

  diffusion(PhiC,dt,Diffv);
}
