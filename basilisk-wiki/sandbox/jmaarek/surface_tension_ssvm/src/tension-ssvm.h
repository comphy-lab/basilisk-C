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

event acceleration(i++){
  // recover the real surface tension if a smooth curvature is used
  for (scalar f in interfaces)
    f.sigma = f.sigma_tmp;
}

#include "tension.h"
#include "curvature-hp.h"

// timestep increase factor
double ratio_dt_st = 1.;

attribute {
  double sigma_tmp;
}

event init (i = 0) {
  for (scalar f in interfaces)
    f.sigma_tmp = f.sigma;
}

/**
## Definition of the potential

We overload the acceleration event to define the potential
$\phi=\sigma\kappa$. Disable this part to the orginal curvature computation. */

event acceleration (i++)
{
  
  /**
  We check for all VOF interfaces for which $\sigma$ is non-zero. */

  for (scalar f in interfaces)
    if (f.sigma) {
      
      /**
      If $\phi$ is already allocated, we add $\sigma\kappa$, otherwise
      we allocate a new field and set it to $\sigma\kappa$. */

      scalar phi = f.phi;
      if (phi.i)
        curvature_hp (f, phi, f.sigma, add = true);
      else {
        phi = new scalar;
        curvature_hp (f, phi, f.sigma, add = false);
        f.phi = phi;
      }

      // disable the accleration event in tension.h
      f.sigma = 0.;
    }
}

event stability(i++) {
  for (scalar f in interfaces)
    f.sigma = f.sigma_tmp/sq(ratio_dt_st);
}

event acceleration(i++) {
  for (scalar f in interfaces)
    f.sigma = f.sigma_tmp;
}

/**
Include any header files here to test the code
*/
// #include "navier-stokes/perfs.h"
#include "viscosity-ssvm.h"
