
extern scalar * tmp_interfaces;

#include "embed_vof.h"

scalar f[], * interfaces = {f};

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;

face vector alphav[];
scalar rhov[];

event defaults (i = 0)
{
  alpha = alphav;
  rho = rhov;
  
  if (mu1 || mu2)
    mu = new face vector;

  display ("draw_vof (c = 'f');");
}

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
#endif

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++)
{

#ifndef sf
for (scalar tmp_c in tmp_interfaces){

  foreach()
    sf[] = (4.*tmp_c[] + 
	    2.*(tmp_c[0,1] + tmp_c[0,-1] + tmp_c[1,0] + tmp_c[-1,0]) +
	    tmp_c[-1,-1] + tmp_c[1,-1] + tmp_c[1,1] + tmp_c[-1,1])/16.;

}
#endif // !sf

#if TREE
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}

event properties (i++)
{
  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    if (mu1 || mu2) {
      face vector muv = mu;
      muv.x[] = fm.x[]*mu(ff);
    }
  }
  
  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}

