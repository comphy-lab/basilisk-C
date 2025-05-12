/**
This is a copy of Lopez' [fracface.h](/sandbox/lopez/src/fracface.h). All credit to him!

Here we extend the function to 3D problems.
*/

#include "geometry.h"
#define VTOL 1.e-6

foreach_dimension()
static double interface_fraction_x (coord m, double alpha, bool right)
{
#if dimension == 2
  alpha += (m.x + m.y)/2;
  coord n = m;
  double xo = (right ? 1. : 0.);
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }
  if (n.y < 1e-4)
    return (n.x*(right ? 1 : -1) < 0. ? 1. : 0.);
  return clamp((alpha - n.x*xo)/n.y, 0., 1.);
#else // dimension == 3
  
  if (fabs(m.y) < 1e-4 && fabs(m.z) < 1e-4)
    return right ? (m.x < 0.) : (m.x > 0.);
  
  double n1, n2;
  double j;
  n1 = m.y/(fabs(m.y) + fabs(m.z));
  n2 = m.z/(fabs(m.y) + fabs(m.z));
  j = right ? 0.5 : -0.5;
  alpha -= j*m.x;
  alpha /= (fabs(m.y) + fabs(m.z));
  return clamp(line_area(n1, n2, alpha), 0., 1.);
#endif  
}

/**
A unique value of the face fraction is calculated from the reconstructed
interfaces at the cells sharing the face by averaging as shown below.
$$
s = \sqrt{vleft \times vright}
$$
*/

void face_fraction (scalar f, face vector s)
{
  boundary({f});

  vector normal_vector[];
  foreach() {
    coord m = mycs (point, f);
    foreach_dimension() 
      normal_vector.x[] = m.x;
  }
  boundary((scalar*){normal_vector});

  foreach_face() {
    if (f[-1] < VTOL || f[] < VTOL) // some cell is empty
      s.x[] = 0.;
    else if (f[-1] > 1.- VTOL && f[] > 1.- VTOL) // both cells are full
      s.x[] = 1.;
    else {
      double vleft = 1., vright = 1.;
      if (f[] < 1. - VTOL) {
        coord m;
	m.x = normal_vector.x[];
	m.y = normal_vector.y[];
	
#if dimension >= 3
	m.z = normal_vector.z[];
#endif
	
        double alpha = plane_alpha (f[], m);

	vleft = interface_fraction_x (m, alpha, false);

      }
      if (f[-1] < 1. - VTOL) {
	coord m;
	m.x = normal_vector.x[-1];
	m.y = normal_vector.y[-1];

#if dimension >= 3
	m.z = normal_vector.z[-1];
#endif
	
        double alpha = plane_alpha (f[-1], m);

	vright = interface_fraction_x (m, alpha, true);
	
      }
      s.x[] = sqrt(vleft*vright);
    }
  }
  boundary((scalar*){s});
}