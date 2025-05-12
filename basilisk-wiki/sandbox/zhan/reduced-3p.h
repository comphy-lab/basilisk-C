/**
# Reduced gravity 

We re-express gravity in [two-phase flows](two-phase.h) as an
[interfacial force](iforce.h) as
$$
-\nabla p + \rho\mathbf{g} = 
-\nabla p' - [\rho]\mathbf{g}\cdot\mathbf{x}\mathbf{n}\delta_s
$$
with $p'= p - \rho\mathbf{g}\cdot\mathbf{x}$ the dynamic pressure and
$\rho\mathbf{g}\cdot\mathbf{x}$ the hydrostatic pressure. The corresponding 
potential is
$$
\phi = [\rho]\mathbf{G}\cdot(\mathbf{x} - \mathbf{Z})
$$
with $\mathbf{G}$ the gravity vector and $\mathbf{Z}$ an optional
reference level. */

coord G = {0.,0.,0.}, Z = {0.,0.,0.};

/**
We need the interfacial force module as well as some
functions to compute the position of the interface. */

#include "iforce.h"
#include "curvature.h"

/**
We overload the acceleration() event to add the contribution of
gravity to the interfacial potential $\phi$.

If $\phi$ is already allocated, we add the contribution of gravity,
otherwise we allocate a new field and set it to the contribution of
gravity. */
  
event acceleration (i++){
  scalar phi1 = f1.phi, phi2 = f2.phi, phi3 = f3.phi;
  coord G1, G2, G3;
  foreach_dimension(){
    G1.x = (- rho1)*G.x;
    G2.x = (- rho2)*G.x;
    G3.x = (- rho3)*G.x;
  }

    
    if (phi1.i)
    position (f1, phi1, G1, Z, add = true);
  else {
    phi1 = new scalar;
    position (f1, phi1, G1, Z, add = false);
    f1.phi = phi1;
  }

  if (phi2.i)
    position (f2, phi2, G2, Z, add = true);
  else {
    phi2 = new scalar;
    position (f2, phi2, G2, Z, add = false);
    f2.phi = phi2;
  }

  if (phi3.i)
    position (f3, phi3, G3, Z, add = true);
  else {
    phi3 = new scalar;
    position (f3, phi3, G3, Z, add = false);
    f3.phi = phi3;
  }
}
