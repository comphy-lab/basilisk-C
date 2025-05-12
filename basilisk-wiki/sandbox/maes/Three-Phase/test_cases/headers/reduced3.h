coord G = {0.,0.,0.}, Z = {0.,0.,0.};


#include "iforce.h"
#include "curvature.h"


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