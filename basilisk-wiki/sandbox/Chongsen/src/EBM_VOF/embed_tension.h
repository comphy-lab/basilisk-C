
#include "embed_iforce.h"
#include "embed_curvature.h"


attribute {
  double sigma;
}

event stability (i++)
{

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin))
    if (fm.x[] > 0.) {
      if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
      if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
      if (Delta < dmin) dmin = Delta;
    }
  double rhom = (1./amin + 1./amax)/2.;


  double sigma = 0.;
  for (scalar c in interfaces)
    sigma += c.sigma;
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(pi*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
}


event acceleration (i++)
{

  for (scalar f in interfaces)
  for (scalar mark in interfaces_mark)
  for (scalar tmp_c in tmp_interfaces){
    if (f.sigma) {
    
      scalar phi = f.phi;
      
      if (phi.i)
  curvature (f, tmp_c, mark, phi, f.sigma, add = true);
      
      else {
  phi = new scalar;
  curvature (f, tmp_c, mark, phi, f.sigma, add = false);
	f.phi = phi;
      }
  
    }
  }

}