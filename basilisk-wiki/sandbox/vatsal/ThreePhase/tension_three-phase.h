/**
# Surface tension

Surface tension can be expressed as the interfacial force density
$$
\phi\nabla f
$$
with $f$ the volume fraction describing the interface and the potential
$$
\phi = \sigma\kappa
$$
with $\sigma$ the (constant) surface tension coefficient and $\kappa$
the interface mean curvature. */

#include "iforce.h"
#include "curvature.h"

/**
The surface tension coefficient is associated to each VOF tracer. */

double sigma12 = 0., sigma23 = 0., sigma13 = 0.;

/**
## Stability condition

The surface tension scheme is time-explicit so the maximum timestep is
the oscillation period of the smallest capillary wave.
$$
T = \sqrt{\frac{\rho_{m}\Delta_{min}^3}{\pi\sigma}}
$$
with $\rho_m=(\rho_1+\rho_2)/2.$ and $\rho_1$, $\rho_2$ the densities
on either side of the interface. */

event stability (i++) {

  /**
  We first compute the minimum and maximum values of $\alpha/f_m =
  1/\rho$, as well as $\Delta_{min}$. */

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin)) {
    if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
    if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
    if (Delta < dmin) dmin = Delta;
  }
  double rhom = (1./amin + 1./amax)/2.;

  /**
  The maximum timestep is set using the sum of surface tension
  coefficients. */

  double sigma = sigma13+sigma23+sigma12;

  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(pi*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
}

/**
## Definition of the potential

We overload the acceleration event to define the potential
$\phi=\sigma\kappa$.<br/>

### This is where we make the change
We first calculate the curvature $\kappa$, instead of $sigma\kappa$. Then we change this field
based on the interface location.
*/

event acceleration (i++)
{
  /**
  We check for all VOF interfaces for which $\sigma$ is non-zero. */
  for (scalar f in interfaces){
    /**
    If $\phi$ is already allocated, an error is thrown. (Fix me!!)
    */
    scalar phi = f.phi;
    if (phi.i){
      fprintf(ferr, "Phi already exists! The current ThreePhase model cannot handle that.\n");
    }
    phi = new scalar;
    curvature (f, phi, 1.0, add = false);
    f.phi = phi;
  }
  /**
  Right now, phi only means $\kappa$. We need to change it to \sigma\kappa based on the interface locations.<br/>
  We refer back to the Schematic:
  <figure>
  <p align="center">
    <img src="https://www.dropbox.com/s/3l0rvxwh6cgrt7g/SchematicBBasilisk.png?dl=1" width="40%">
    <figcaption><p align="center">Schematic: Details of the VOF fields and interfaces.</figcaption>
  </figure>
  * If both phi1[] and phi2[] are defined, then it is the interface between phase 2 and 3.
  * If only phi1[] is defined (but phi2[]) is not, then it is the interface between 1 and 3.
  * If only phi2[] is defined (but phi1[]) is not, then it is the interface between 1 and 2.
  We modify $phi$ accordingly.
  */
  scalar phi1 = f1.phi, phi2 = f2.phi;
  foreach(){ // To do: Use a better if-else look. Also check if for (scalar f in interfaces) can be used somehow.
    if(phi1[] < nodata && phi2[] < nodata){
      phi1[] *= sigma23/2.;
      phi2[] *= sigma23/2.;
    } else {
      if(phi1[] < nodata){
        phi1[] *= sigma13;
      } else {
        if (phi2[] < nodata){
          phi2[] *= sigma12;
        }
      }
    }
  }
}
