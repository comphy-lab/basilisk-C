/**
#Surface tension
This file aims to add surface tension to the multilayer solver.
The surface tension field $\sigma$ is introduced. It must be initialized to add surface tension. 
It can be a function, especially of liquid composition or solute concentration. */
#include "curvature.h"
static bool Laplace = true;
static bool nh_tension = true;
static bool corr_dux = false;
static bool corr_dwx = false;
scalar sigma[];

/**
## Stability condition
The surface tension scheme is time-explicit so the maximum timestep is the oscillation period of 
the smallest capillary wave. The maximum timestep is set using the standard (hydrostatic) CFL condition
$$
\Delta t < \text{CFL}\frac{\Delta}{|\mathbf{u}|+\sqrt{gH + \frac{\sigma\,H}{\rho\,\Delta}}}
$$
or the modified (non-hydrostatic) condition([Popinet,
2020](/Bibliography#popinet2020)), corrected with capilary waves
$$
\Delta t < 
\text{CFL}\frac{\Delta}{\left|\vec{u}\right| + \sqrt{\pp{g\,\Delta + \dfrac{\sigma}{\rho\,\Delta}}\tanh{\pp{\dfrac{H}{\Delta}}}}}
$$
*/
event stability (i++)
{
  double rhom = 1.;

  foreach_face (reduction (min:dtmax)) {
    double Hf = 0.;
    foreach_layer()
      Hf += hf.x[];
    if (Hf > dry) {
      Hf /= fm.x[];
      if (sigma[] > 0.) {
      double cp = sqrt((fabs(G)*Delta + pi*sigma[]/(rhom*Delta))*tanh(Hf/Delta));            
        foreach_layer() {
          double c = fabs(hu.x[]/hf.x[])/CFL + cp/CFL_H;
          if (c > 0){
            double dt = min(cm[], cm[-1])*Delta/(c*fm.x[]);
         	  if (dt < dtmax)
       	      dtmax = dt;
          }
      	}
      }
    }
  }
}

/**
## Marangoni flows */
event viscous_term (i++) 
{
  if (nu > 0) {
    vector du[];
    foreach()
      foreach_dimension() {
        du.x[] = - corr_dwx*(w[1,0,nl-1]-w[-1,0,nl-1])/(2*Delta) + (sigma[1] - sigma[-1])/(2*Delta*hm(point, true)*nu);
        }
    dut = du;
  }
}

/**
## Laplace pressure
The calculation of the acceleration is done by this event, overloaded from 
[its definition](layered/hydro.h) in the multilayer Navier--Stokes solver. 
*/
event pressure (i++) 
{
  if (nh_tension == true) {
    foreach()
      phi_s[] = - Laplace * sigma[] * kappa_c(point) - corr_dux * 2.*nu*(u.x[1,0,nl-1]-u.x[-1,0,nl-1])/(2*Delta);
    boundary({phi_s});
  }
  else {
    scalar phi_t[];
    foreach()
      phi_t[] = - Laplace * sigma[] * kappa_c(point);
    boundary({phi_t});
    
    foreach_face()
      foreach_layer()
        if (hf.x[] > dry) {
        	double ax = fm.x[]*(phi_t[] -	phi_t[-1])/Delta;
        	hu.x[] -= hf.x[]*dt*ax;
        	ha.x[] -= hf.x[]*ax;
        }
    boundary ((scalar *){hu});
    boundary_flux ({ha}); 
    }     
}